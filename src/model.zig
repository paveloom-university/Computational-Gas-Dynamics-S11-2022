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
                .writer = writer,
            };
        }
        /// Compute the flow to/from the neighbour cell
        inline fn flow(
            self: *Self,
            comptime dir: Dir,
            vel: anytype,
            ro: anytype,
            index: usize,
        ) F {
            if (self.grid.nearEdge(dir, index)) {
                return 0;
            } else {
                const border_value = self.grid.border(dir, vel, index);
                return self.h * self.tau * border_value * if (border_value > 0)
                    ro[index]
                else
                    self.grid.neighbour(dir, ro, index);
            }
        }
        /// Compute an update to the velocities and energy on the final stage
        inline fn final_update(
            self: *Self,
            array: anytype,
            index: usize,
            ro: F,
            ro_prev: F,
            f_top: F,
            f_bottom: F,
            f_left: F,
            f_right: F,
        ) F {
            return ro_prev / ro * array[index] +
                (self.grid.neighbour(.top, array, index) * f_top +
                self.grid.neighbour(.bottom, array, index) * f_bottom +
                self.grid.neighbour(.left, array, index) * f_left +
                self.grid.neighbour(.right, array, index) * f_right) /
                ro / self.h / self.h;
        }
        /// Compute the evolution of the system for a single time step
        fn step(self: *Self, slice: *Cells(F).Slice) void {
            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Stage 1");
            defer zone.end();

            // Unpack the struct data
            const tau = self.tau;
            const h = self.h;
            const m = self.grid.n * self.grid.n;
            // Get the fields of interest from the slice
            const u = slice.items(.u);
            const v = slice.items(.v);
            const ro = slice.items(.ro);
            const e = slice.items(.ro);
            const p = slice.items(.p);
            const u_aux = slice.items(.u_aux);
            const v_aux = slice.items(.v_aux);
            const e_aux = slice.items(.e_aux);
            // For each cell in the grid
            var index: usize = 0;
            while (index < m) : (index += 1) {
                // Execute the Euler stage
                //
                // At this stage, only the quantities related to the cell as
                // a whole change, and the liquid is assumed to be retarded.

                // Compute the auxiliary values at the borders of the cells
                const p_top = self.grid.border(.top, p, index);
                const p_bottom = self.grid.border(.bottom, p, index);
                const p_left = self.grid.border(.left, p, index);
                const p_right = self.grid.border(.right, p, index);
                const v_top = self.grid.border(.top, v, index);
                const v_bottom = self.grid.border(.bottom, v, index);
                const u_left = self.grid.border(.left, u, index);
                const u_right = self.grid.border(.right, u, index);
                // Compute the auxiliary values inside the cells
                const k = tau / h / ro[index];
                u_aux[index] = u[index] - k * (p_left - p_right);
                v_aux[index] = v[index] - k * (p_top - p_bottom);
                e_aux[index] = e[index] - k * (p_top * v_top - p_bottom * v_bottom + p_left * u_left - p_right * u_right);
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
                const f_top = self.flow(.top, v_aux, ro, index);
                const f_bottom = self.flow(.bottom, v_aux, ro, index);
                const f_left = self.flow(.left, u_aux, ro, index);
                const f_right = self.flow(.right, u_aux, ro, index);

                // Execute the final stage
                //
                // Here, the redistribution of mass, momentum, and energy
                // over space occurs and the final fields of the Euler
                // flow parameters are determined

                const ro_prev = ro[index];
                ro[index] += (f_top + f_bottom + f_left + f_right) / h / h;
                u[index] = self.final_update(u_aux, index, ro[index], ro_prev, f_top, f_bottom, f_left, f_right);
                v[index] = self.final_update(v_aux, index, ro[index], ro_prev, f_top, f_bottom, f_left, f_right);
                e[index] = self.final_update(e_aux, index, ro[index], ro_prev, f_top, f_bottom, f_left, f_right);
            }
        }
        /// Compute the evolution of the system for a specific amount of time steps
        pub fn compute(self: *Self, s: usize) !void {
            // Compute pointers to the start of each field of the array of cells
            var slice = self.grid.cells.slice();
            // Write the header
            try self.writer.writeHeader(s, self.grid.n);
            // For each time step
            var i: usize = 0;
            while (i < s) : (i += 1) {
                // Perform a computation step
                self.step(&slice);
                // Write the current grid to the output file
                try self.writer.writeGrid(i, &self.grid);
            }
        }
    };
}
