const std = @import("std");

const tracy = @import("tracy");

const grid = @import("grid.zig");

const Cells = grid.Cells;
const Dir = grid.Dir;
const Grid = grid.Grid;
const Path = @import("lib.zig").Path;
const Writer = @import("writer.zig").Writer;

const tracker = 200;

/// Model
///
/// Note the representation of cells in memory doesn't change
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
        tau: F,
        h: F,
        grid: Grid(F),
        eqs: fn (u: F, v: F, ro: F, e: F) F,
        /// The writer is initialized from the path
        writer: Writer(F),
        /// Initialize the model
        pub fn init(s: struct {
            /// Time step
            tau: F,
            /// Grid step
            h: F,
            /// The Euler grid
            grid: Grid(F),
            /// Equation of state
            eqs: fn (
                /// Velocity component along the X axis
                u: F,
                /// Velocity component along the X axis
                v: F,
                /// Density
                ro: F,
                /// Specific energy
                e: F,
            ) F,
            /// Relative paths to outputs files
            ///
            /// The results will be written to these each time step.
            /// Note that only extensions `.csv` (for text files)
            /// and `.bin` (for binary files) are allowed, and
            /// there must be at least one file
            path: Path,
        }) !Self {
            // Prepare the writer
            const writer = try Writer(F).from(s.path);
            // Initialize the model
            return Self{
                .tau = s.tau,
                .h = s.h,
                .grid = s.grid,
                .eqs = s.eqs,
                .writer = writer,
            };
        }
        /// A flow to/from the neighbour cell
        const Flow = struct {
            value: F,
            /// 1.0 if the gas flows into the cell, 0.0 otherwise
            into: F,
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
            const aux_vel_average = (vel_aux[index] + self.grid.neighbour(dir, true, vel_aux, index)) / 2;
            if (index == tracker + self.grid.n or index == tracker) {
                std.debug.print("aux_vel_average: {} {}\n", .{ index, aux_vel_average });
            }
            // const into: F = if (aux_vel_average > 0) 0.0 else 1.0;
            const into: F = if (aux_vel_average > 0) switch (dir) {
                .right, .top => 0.0,
                .left, .bottom => 1.0,
            } else switch (dir) {
                .right, .top => 1.0,
                .left, .bottom => 0.0,
            };
            // const value = aux_vel_average * if (aux_vel_average > 0)
            //     ro_aux[index]
            // else
            //     self.grid.neighbour(dir, false, ro_aux, index);
            const value = aux_vel_average * ((1 - into) * ro_aux[index] +
                into * self.grid.neighbour(dir, false, ro_aux, index));
            // const value = aux_vel_average * if (aux_vel_average > 0) switch (dir) {
            //     .right, .top => ro_aux[index],
            //     .left, .bottom => self.grid.neighbour(dir, false, ro_aux, index),
            // } else switch (dir) {
            //     .right, .top => -self.grid.neighbour(dir, false, ro_aux, index),
            //     .left, .bottom => -ro_aux[index],
            // };
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
            const first = ro_prev / ro * array[index];
            // const second = (flows.left.value * (1 - flows.left.into) * self.grid.neighbour(.left, normal, array, index) -
            //     flows.right.value * (1 - flows.right.into) * self.grid.neighbour(.right, normal, array, index) +
            //     flows.bottom.value * (1 - flows.bottom.into) * self.grid.neighbour(.bottom, normal, array, index) -
            //     flows.top.value * (1 - flows.top.into) * self.grid.neighbour(.top, normal, array, index) +
            //     (flows.left.value * flows.left.into -
            //     flows.right.value * flows.right.into +
            //     flows.bottom.value * flows.bottom.into -
            //     flows.top.value * flows.top.into) * array[index]) *
            //     self.tau / self.h / ro;
            const second = (flows.left.value * flows.left.into * self.grid.neighbour(.left, normal, array, index) -
                flows.right.value * flows.right.into * self.grid.neighbour(.right, normal, array, index) +
                flows.bottom.value * flows.bottom.into * self.grid.neighbour(.bottom, normal, array, index) -
                flows.top.value * flows.top.into * self.grid.neighbour(.top, normal, array, index) +
                (flows.left.value * (1 - flows.left.into) -
                flows.right.value * (1 - flows.right.into) +
                flows.bottom.value * (1 - flows.bottom.into) -
                flows.top.value * (1 - flows.top.into)) * array[index]) *
                self.tau / self.h / ro;
            if (index == tracker) {
                std.debug.print("final_update: {} {}\n", .{ first, second });
            }
            return first + second;
        }
        /// Compute the evolution of the system for a single time step
        fn step(self: *Self, slice: *Cells(F).Slice) void {
            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Step");
            defer zone.end();

            // Unpack the struct data
            const tau = self.tau;
            const h = self.h;
            const m = self.grid.n * self.grid.n;
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
            var index: usize = 0;
            while (index < m) : (index += 1) {
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
                const k = tau / h / ro[index];
                u_aux[index] = u[index] - k * (p_right - p_left);
                v_aux[index] = v[index] - k * (p_top - p_bottom);
                e_aux[index] = e[index] - k * (p_top * v_top - p_bottom * v_bottom + p_right * u_right - p_left * u_left);
                if (index == tracker) {
                    std.debug.print("{}\n", .{index});
                    std.debug.print("u_aux: {}, v_aux: {} e_aux: {}\n", .{ u_aux[index], v_aux[index], e_aux[index] });
                    std.debug.print("k: {}\n", .{k});
                    std.debug.print("p_top: {}, p_bottom: {}, p_left: {}, p_right: {}\n", .{ p_top, p_bottom, p_left, p_right });
                    std.debug.print("v_top: {}, v_bottom: {}, u_left: {}, u_right: {}\n", .{ v_top, v_bottom, u_left, u_right });
                    // std.debug.print("ro: {}, u: {}, v: {}, e: {}\n", .{ ro[index], u[index], v[index], e[index] });
                }
            }
            // For each cell in the grid
            index = 0;
            while (index < m) : (index += 1) {
                // Copy the density matrix (since we compute
                // the density during the time step)
                ro_aux[index] = ro[index];
            }
            // For each cell in the grid
            index = 0;
            while (index < m) : (index += 1) {
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
                ro[index] += tau / h * (flows.left.value - flows.right.value + flows.bottom.value - flows.top.value);
                // ro[index] += 2 * tau / h *
                //     (flows.bottom.value * (flows.bottom.into - 0.5) +
                //     flows.left.value * (flows.left.into - 0.5) +
                //     flows.top.value * (flows.top.into - 0.5) +
                //     flows.right.value * (flows.right.into - 0.5));
                u[index] = self.final_update(true, u_aux, index, ro[index], ro_prev, flows);
                v[index] = self.final_update(true, v_aux, index, ro[index], ro_prev, flows);
                e[index] = self.final_update(false, e_aux, index, ro[index], ro_prev, flows);
                p[index] = self.eqs(u[index], v[index], ro[index], e[index]);
                if (index == tracker + self.grid.n) {
                    std.debug.print("{} flows: {}\n\n", .{ index, flows });
                }
                if (index == tracker) {
                    std.debug.print("flows: {}\n", .{flows});
                    // std.debug.print("u_aux: {}, v_aux: {} e_aux: {}\n", .{ u_aux[index], v_aux[index], e_aux[index] });
                    // std.debug.print("{} {} {} {} {} {} {}\n", .{ p_top, p_bottom, p_left, p_right, p[index], p[index - self.grid.n], k });
                    std.debug.print("ro: {}, u: {}, v: {}, e: {}\n", .{ ro[index], u[index], v[index], e[index] });
                }
            }
        }
        /// Compute the evolution of the system for a specific amount of time steps
        pub fn compute(self: *Self, s: usize) !void {
            // Compute pointers to the start of each field of the array of cells
            var slice = self.grid.cells.slice();
            // Write the header
            try self.writer.writeHeader(s, self.grid.n);
            // Write the initial grid to the output file
            try self.writer.writeGrid(&slice);
            // For each time step
            var i: usize = 0;
            while (i < s) : (i += 1) {
                // Perform a computation step
                self.step(&slice);
                // Write the current grid to the output file
                try self.writer.writeGrid(&slice);
            }
        }
    };
}
