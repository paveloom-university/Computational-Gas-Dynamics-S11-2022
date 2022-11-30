//! An implementation of the large-particle method for the numerical
//! solution of hydrodynamic equations. Specifically, focusing on
//! solving a convection equation

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
        const Self = @This();
        /// Velocity component along the X axis
        u: F,
        /// Velocity component along the X axis
        v: F,
        /// Density
        ro: F,
        /// Specific energy (energy per unit mass)
        e: F,
        /// Pressure (as a function of density and specific inner energy)
        p: F,
        /// Auxiliary velocity component along the X axis
        u_aux: F = 0,
        /// Auxiliary velocity component along the X axis
        v_aux: F = 0,
        /// Auxiliary specific energy (energy per unit mass)
        e_aux: F = 0,
        /// Create a data record from a cell (as bytes)
        fn record(self: *const Self) [32]u8 {
            return std.mem.toBytes([_]f64{ self.u, self.v, self.ro, self.e });
        }
    };
}

// Cells in the Euler grid
pub fn Cells(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return std.MultiArrayList(Cell(F));
}

/// Euler grid
///
/// The returned struct emulates an N x N matrix of cells
pub fn Grid(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return struct {
        /// Size of the grid
        n: usize,
        /// The underlying data
        cells: Cells(F),
    };
}

/// Path
const Path = []const u8;

/// Output file writer
fn Writer(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return struct {
        const Self = @This();
        const FileWriter = std.fs.File.Writer;
        const BufferedWriter = std.io.BufferedWriter(4096, FileWriter);
        writer: BufferedWriter,
        inline fn writeUsize(self: *Self, buf: *[@sizeOf(usize)]u8, value: usize) !void {
            std.mem.writeIntNative(usize, buf, value);
            _ = try self.writer.write(buf);
        }
        /// Write the header
        pub fn writeHeader(self: *Self, s: usize, n: usize) !void {
            // Prepare a buffer for the bytes
            var buf: [@sizeOf(usize)]u8 = undefined;
            // Write the number of time steps, plus one
            try self.writeUsize(&buf, s + 1);
            // Write the size of the grid
            try self.writeUsize(&buf, n);
        }
        /// Write the grid
        pub fn writeGrid(self: *Self, i: usize, grid: *const Grid(F)) !void {
            // If it's the first step
            if (i == 0) {}
            // For each cell
            var j: usize = 0;
            while (j < grid.n * grid.n) : (j += 1) {
                const cell = grid.cells.get(j);
                // Create a data record
                const record = cell.record();
                // Write the data
                _ = try self.writer.write(&record);
            }
        }
        /// Create a writer from the path
        pub fn from(path: Path) !Self {
            // Create a file
            const file = if (std.fs.path.isAbsolute(path))
                try std.fs.createFileAbsolute(path, .{})
            else
                try std.fs.cwd().createFile(path, .{});
            // Return the writer
            return Self{ .writer = std.io.bufferedWriter(file.writer()) };
        }
    };
}

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
