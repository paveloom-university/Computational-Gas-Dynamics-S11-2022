const tracy = @import("tracy");

const grid = @import("grid.zig");

const Cells = grid.Cells;
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
        /// Compute a value at the top border of the cell
        inline fn top(array: anytype, n: usize, index: usize) F {
            return if (index / n == 0)
                array[index]
            else
                (array[index] + array[index - n]) / 2;
        }
        /// Compute a value at the bottom border of the cell
        inline fn bottom(array: anytype, n: usize, index: usize) F {
            return if (index / n == n - 1)
                array[index]
            else
                (array[index] + array[index + n]) / 2;
        }
        /// Compute a value at the left border of the cell
        inline fn left(array: anytype, n: usize, index: usize) F {
            return if (index % n == 0)
                array[index]
            else
                (array[index] + array[index - 1]) / 2;
        }
        /// Compute a value at the right border of the cell
        inline fn right(array: anytype, n: usize, index: usize) F {
            return if (index % n == n - 1)
                array[index]
            else
                (array[index] + array[index + 1]) / 2;
        }
        /// Execute stage 1: the Euler step
        ///
        /// At this stage, only the quantities related to the cell as
        /// a whole change, and the liquid is assumed to be retarded.
        fn stage_1(self: *Self, slice: *Cells(F).Slice) void {
            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Stage 1");
            defer zone.end();
            // Unpack the struct data
            const tau = self.tau;
            const h = self.h;
            const n = self.grid.n;
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
            var i: usize = 0;
            var j: usize = 0;
            while (j < n) : (j += 1) {
                while (i < n) : (i += 1) {
                    const index = j * n + i;
                    // Compute the auxiliary values on the borders of the cells
                    const p_top = top(p, n, index);
                    const p_bottom = bottom(p, n, index);
                    const p_left = left(p, n, index);
                    const p_right = right(p, n, index);
                    const v_top = top(v, n, index);
                    const v_bottom = bottom(v, n, index);
                    const u_left = left(u, n, index);
                    const u_right = right(u, n, index);
                    // Compute the auxiliary values inside the cells
                    const k = tau / h / ro[index];
                    u_aux[index] = u[index] - k * (p_left - p_right);
                    v_aux[index] = v[index] - k * (p_top - p_bottom);
                    e_aux[index] = e[index] - k *
                        (p_top * v_top - p_bottom * v_bottom +
                        p_left * u_left - p_right * u_right);
                }
            }
        }
        /// Compute the evolution of the system for a single time step
        fn step(self: *Self, slice: *Cells(F).Slice) void {
            // Execute all stages
            self.stage_1(slice);
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
