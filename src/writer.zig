const std = @import("std");

const Grid = @import("grid.zig").Grid;
const Path = @import("lib.zig").Path;

/// Output file writer
pub fn Writer(
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
