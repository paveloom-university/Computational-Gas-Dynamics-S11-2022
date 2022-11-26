//! An implementation of the large-particle method for the numerical
//! solution of hydrodynamic equations. Specifically, focusing on
//! the evolution of the motion of a calorically perfect gas.

const std = @import("std");

const tracy = @import("tracy");

/// A cell in the Euler grid
pub fn Cell(
    /// Type of a floating-point number
    comptime F: type,
) type {
    // Check if the provided type isn't a floating-point number type
    if (@typeInfo(F) != .Float) {
        @compileError("This type isn't allowed for floats");
    }
    return struct {
        /// Velocity component along the X axis
        u: F = 0,
        /// Velocity component along the X axis
        v: F = 0,
        /// Density
        ro: F = 0,
        /// Specific energy (energy per unit mass)
        e: F = 0,
        /// Pressure (as a function of density and specific inner energy)
        p: F = 0,
        /// Auxiliary velocity component along the X axis
        u_aux: F = 0,
        /// Auxiliary velocity component along the X axis
        v_aux: F = 0,
        /// Auxiliary specific energy (energy per unit mass)
        e_aux: F = 0,
    };
}

/// Euler grid
///
/// The returned struct emulates an N x N matrix of cells
pub fn Grid(
    /// Type of a floating-point number
    comptime F: type,
) type {
    // Check if the provided type isn't a floating-point number type
    if (@typeInfo(F) != .Float) {
        @compileError("This type isn't allowed for floats");
    }
    return struct {
        const Self = @This();
        /// Type of the underlying data
        const Cells = std.MultiArrayList(Cell(F));
        /// Allocator
        allocator: std.mem.Allocator,
        /// Time step
        tau: F,
        /// Grid step
        h: F,
        /// Size of the grid
        n: usize,
        /// The underlying data
        cells: Cells,
        /// Initialize a grid with a specific size
        pub fn init(
            args: struct {
                /// Allocator
                allocator: std.mem.Allocator,
                /// Time step
                tau: F,
                /// Grid step
                h: F,
                /// Size of the grid
                n: usize,
            },
        ) !Self {
            // Initialize the underlying data
            var cells = Cells{};
            // Compute the number of cells
            const m = args.n * args.n;
            // Prepare space for storing N x N cells
            try cells.ensureTotalCapacity(args.allocator, m);
            // Initialize each cell
            var i: usize = 0;
            while (i < m) : (i += 1) {
                cells.appendAssumeCapacity(Cell(F){});
            }
            // Return the grid
            return Self{
                .allocator = args.allocator,
                .tau = args.tau,
                .h = args.h,
                .n = args.n,
                .cells = cells,
            };
        }
        /// Deinitialize the grid
        pub fn deinit(self: *Self) void {
            self.cells.deinit(self.allocator);
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
        fn stage_1(self: *Self, slice: *Cells.Slice) void {
            // Mark the Tracy zone
            const zone = tracy.ZoneN(@src(), "Stage 1");
            defer zone.end();
            // Unpack the struct data
            const tau = self.tau;
            const h = self.h;
            const n = self.n;
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
        fn step(self: *Self, slice: *Cells.Slice) void {
            // Execute all stages
            self.stage_1(slice);
        }
        /// Compute the evolution of the system for a specific amount of time steps
        pub fn compute(self: *Self, n: usize) void {
            // Compute pointers to the start of each field of the array of cells
            var slice = self.cells.slice();
            // For each time step
            var i: usize = 0;
            while (i < n) : (i += 1) {
                // Perform a computation step
                self.step(&slice);
            }
        }
    };
}
