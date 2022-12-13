const std = @import("std");

const tracy = @import("tracy");

const grid = @import("grid.zig");

const Cells = grid.Cells;
const Dir = grid.Dir;
const Grid = grid.Grid;
const Path = @import("lib.zig").Path;
const Writer = @import("writer.zig").Writer;

/// Model
///
/// Note that the representation of the cells in memory doesn't change
/// with the grid step. Only the latter is model-specific.
pub fn Model(
    /// Type of a floating-point number
    comptime F: type,
) type {
    // Check if the provided type isn't a floating-point number type
    if (@typeInfo(F) != .Float) {
        @compileError("This type isn't allowed for floats");
    }
    return struct {
        const Self = @This();
        /// A collection of async frames. These are used for multithreading
        const Promises = struct {
            step_1: []@Frame(step1),
            step_2: []@Frame(step2),
        };
        allocator: std.mem.Allocator,
        tau: F,
        h: F,
        grid: Grid(F),
        phi: F,
        eqs: fn (u: F, v: F, ro: F, e: F) F,
        /// The writer is initialized from the path
        writer: Writer(F),
        threads: usize,
        promises: Promises = undefined,
        /// Initialize the model
        pub fn init(s: struct {
            /// Allocator (you're expected to allocate `Grid.Cells` yourself,
            /// but the allocator is also used for allocating promises)
            allocator: std.mem.Allocator,
            /// Time step [s]
            tau: F,
            /// Grid step [m]
            h: F,
            /// The Euler grid
            grid: Grid(F),
            /// The value of the gravity potential [m^2/s^2]
            ///
            /// We expect it to be a constant for all cells
            phi: F,
            /// Equation of state
            eqs: fn (
                /// Velocity component along the X axis
                u: F,
                /// Velocity component along the Y axis
                v: F,
                /// Density
                ro: F,
                /// Specific energy
                e: F,
            ) F,
            /// Relative path to the output file
            path: Path,
            /// Number of CPU threads to use
            threads: usize,
        }) !Self {
            // Prepare the writer
            const writer = try Writer(F).from(s.path);
            // Initialize the model
            var self = Self{
                .allocator = s.allocator,
                .tau = s.tau,
                .h = s.h,
                .grid = s.grid,
                .eqs = s.eqs,
                .phi = s.phi,
                .writer = writer,
                .threads = s.threads,
            };
            // If there are several threads
            if (self.threads > 1) {
                // Prepare arrays for promises
                self.promises = Promises{
                    .step_1 = try s.allocator.alloc(@Frame(step1), s.threads),
                    .step_2 = try s.allocator.alloc(@Frame(step2), s.threads),
                };
            }
            return self;
        }
        /// Deinitialize the model (this includes
        /// the `Grid.Cells` and the promises)
        pub fn deinit(self: *Self) void {
            self.grid.cells.deinit(self.allocator);
            if (self.threads > 1) {
                self.allocator.free(self.promises.step_1);
                self.allocator.free(self.promises.step_2);
            }
            self.* = undefined;
        }
        /// A flow to/from the neighbour cell
        const Flow = struct {
            value: F,
            /// `true` if the gas flows into the cell, `false` otherwise
            into: bool,
        };
        /// Flows to/from the neighbour cells
        const Flows = struct {
            top: Flow,
            bottom: Flow,
            left: Flow,
            right: Flow,
        };
        /// Compute the flow to/from the neighbour cell
        inline fn flow(
            self: *Self,
            comptime dir: Dir,
            vel_aux: []F,
            ro_aux: []F,
            index: usize,
        ) Flow {
            // Compute the halved sum of the velocities
            const aux_vel_average = (vel_aux[index] + self.grid.neighbour(dir, true, vel_aux, index)) / 2;
            // Since we have axes left-to-right, bottom-to-top,
            // having positive sum of velocities here means
            // that the gas flows out of the cell if the direction
            // is right or top, and into the cell otherwise.
            const into = if (aux_vel_average > 0) switch (dir) {
                .right, .top => false,
                .left, .bottom => true,
            } else switch (dir) {
                .right, .top => true,
                .left, .bottom => false,
            };
            // Whether gas flows into the cell or out of it defines whose density we take
            const value =
                (if (into) self.grid.neighbour(dir, false, ro_aux, index) else ro_aux[index]) *
                aux_vel_average * self.tau * self.h;
            // Return the flow
            return Flow{
                .value = value,
                .into = into,
            };
        }
        /// Compute an update to the velocities and energy on the final stage
        inline fn final_update(
            self: *Self,
            comptime normal: bool,
            array: []F,
            index: usize,
            ro: F,
            ro_prev: F,
            flows: Flows,
        ) F {
            return ro_prev / ro * array[index] +
                // Note that we choose between the neighbour and the current cell here, too
                (flows.left.value * (if (flows.left.into) self.grid.neighbour(.left, normal, array, index) else array[index]) -
                flows.right.value * (if (flows.right.into) self.grid.neighbour(.right, normal, array, index) else array[index]) +
                flows.bottom.value * (if (flows.bottom.into) self.grid.neighbour(.bottom, normal, array, index) else array[index]) -
                flows.top.value * (if (flows.top.into) self.grid.neighbour(.top, normal, array, index) else array[index])) /
                self.h / self.h / ro;
        }
        /// Compute the first part of the time step
        fn step1(self: *Self, slice: *Cells(F).Slice, start: usize, end: usize) void {
            // If there are several threads
            if (self.threads > 1) {
                // Inform the event loop we're CPU bound.
                // This effectively puts a worker on every logical core.
                std.event.Loop.startCpuBoundOperation();
            }

            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Step 1");
            defer zone.end();

            // Unpack the model parameters
            const tau = self.tau;
            const h = self.h;
            const phi = self.phi;
            // Get the fields of interest from the slice
            const u = slice.items(.u);
            const v = slice.items(.v);
            const ro = slice.items(.ro);
            const e = slice.items(.e);
            const p = slice.items(.p);
            const u_aux = slice.items(.u_aux);
            const v_aux = slice.items(.v_aux);
            const e_aux = slice.items(.e_aux);
            // For each cell in the grid
            var index: usize = start;
            while (index < end) : (index += 1) {
                // Execute the Euler stage
                //
                // At this stage, only the quantities related to the cell as
                // a whole change, and the liquid is assumed to be retarded.

                // Compute the auxiliary values at the borders of the cells
                const p_top = self.grid.border(.top, false, p, index);
                const p_bottom = self.grid.border(.bottom, false, p, index);
                const p_left = self.grid.border(.left, false, p, index);
                const p_right = self.grid.border(.right, false, p, index);
                const v_top = self.grid.border(.top, true, v, index);
                const v_bottom = self.grid.border(.bottom, true, v, index);
                const u_left = self.grid.border(.left, true, u, index);
                const u_right = self.grid.border(.right, true, u, index);
                // Compute the auxiliary values inside the cells
                const d = tau / h;
                const k = d / ro[index];
                u_aux[index] = u[index] - k * (p_right - p_left);
                // We expect this extra term to to act as a gravity potential
                v_aux[index] = v[index] - k * (p_top - p_bottom) - d * phi;
                e_aux[index] = e[index] - k * (p_top * v_top - p_bottom * v_bottom + p_right * u_right - p_left * u_left);
            }
        }
        /// Compute the second part of the time step
        fn step2(self: *Self, slice: *Cells(F).Slice, start: usize, end: usize) void {
            // If there are several threads
            if (self.threads > 1) {
                // Inform the event loop we're CPU bound.
                // This effectively puts a worker on every logical core.
                std.event.Loop.startCpuBoundOperation();
            }

            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Step 2");
            defer zone.end();

            // Unpack the model parameters
            const h = self.h;
            // Get the fields of interest from the slice
            const u = slice.items(.u);
            const v = slice.items(.v);
            const ro = slice.items(.ro);
            const e = slice.items(.e);
            const p = slice.items(.p);
            const u_aux = slice.items(.u_aux);
            const v_aux = slice.items(.v_aux);
            const ro_aux = slice.items(.ro_aux);
            const e_aux = slice.items(.e_aux);
            // For each cell in the grid
            var index: usize = start;
            while (index < end) : (index += 1) {
                // Execute the Lagrangian stage
                //
                // At this stage, transfer effects are calculated that take
                // into account the exchange between cells when they are
                // rebuilt to the previous Euler grid

                // Note that we expect no incoming flow outside of the grid
                const flows = Flows{
                    .top = self.flow(.top, v_aux, ro_aux, index),
                    .bottom = self.flow(.bottom, v_aux, ro_aux, index),
                    .left = self.flow(.left, u_aux, ro_aux, index),
                    .right = self.flow(.right, u_aux, ro_aux, index),
                };

                // Execute the final stage
                //
                // Here, the redistribution of mass, momentum, and energy
                // over space occurs and the final fields of the Euler
                // flow parameters are determined

                const ro_prev = ro[index];
                ro[index] += (flows.left.value - flows.right.value + flows.bottom.value - flows.top.value) / h / h;
                u[index] = self.final_update(true, u_aux, index, ro[index], ro_prev, flows);
                v[index] = self.final_update(true, v_aux, index, ro[index], ro_prev, flows);
                e[index] = self.final_update(false, e_aux, index, ro[index], ro_prev, flows);
                p[index] = self.eqs(u[index], v[index], ro[index], e[index]);
            }
        }
        /// Compute the evolution of the system for a single time step
        fn step(self: *Self, slice: *Cells(F).Slice) !void {
            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Step");
            defer zone.end();

            // Get the fields of interest from the slice
            const ro = slice.items(.ro);
            const ro_aux = slice.items(.ro_aux);
            // Divide cells between the threads
            const m = self.grid.n * self.grid.n;
            const h = m / self.threads;
            const r = m % self.threads;
            // If there are several threads
            if (self.threads > 1) {
                // For each thread
                var thread: usize = 0;
                while (thread < self.threads) : (thread += 1) {
                    // The first thread here may take a bit more work than others
                    const start = if (thread == 0) 0 else r + h * thread;
                    const end = r + h * (thread + 1);
                    // Execute the first part of the method
                    // (includes the Euler stage)
                    self.promises.step_1[thread] = async self.step1(slice, start, end);
                }
                // Synchronize the threads
                for (self.promises.step_1) |*future| {
                    await future;
                }
            } else {
                // Otherwise, just do a simple call
                self.step1(slice, 0, m);
            }
            // Copy the density matrix (since we
            // compute the density in the next part)
            std.mem.copy(F, ro_aux, ro);
            // If there are several threads
            if (self.threads > 1) {
                // For each thread
                var thread: usize = 0;
                while (thread < self.threads) : (thread += 1) {
                    // The first thread here may take a bit more work than others
                    const start = if (thread == 0) 0 else r + h * thread;
                    const end = r + h * (thread + 1);
                    // Execute the second part of the method
                    // (includes the Lagrangian stage and the final stage)
                    self.promises.step_2[thread] = async self.step2(slice, start, end);
                }
                // Synchronize the threads
                for (self.promises.step_2) |*future| {
                    await future;
                }
            } else {
                // Otherwise, just do a simple call
                self.step2(slice, 0, m);
            }
        }
        /// Compute the evolution of the system for a specific amount
        /// of time steps `s`, saving the results every `d` frames
        pub fn compute(self: *Self, s: usize, d: usize) !void {
            // Compute pointers to the start of each field of the array of cells
            var slice = self.grid.cells.slice();
            // Write the header
            try self.writer.writeHeader(s, d, self.grid.n);
            // Write the initial grid to the output file
            try self.writer.writeGrid(&slice);
            // For each time step
            var i: usize = 1;
            while (i <= s) : (i += 1) {
                // Perform a computation step
                try self.step(&slice);
                // Save the results every `d` frames, but
                // make sure the last one is saved, too
                if (i % d == 0 or i == s) {
                    // Write the current grid to the output file
                    try self.writer.writeGrid(&slice);
                }
            }
        }
    };
}
