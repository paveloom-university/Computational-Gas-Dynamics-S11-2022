const std = @import("std");

const grid = @import("grid.zig");

const Cells = grid.Cells;
const Grid = grid.Grid;
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
        fn writeUsize(self: *Self, buf: *[@sizeOf(usize)]u8, value: usize) !void {
            std.mem.writeIntNative(usize, buf, value);
            _ = try self.writer.write(buf);
        }
        /// Write the header
        pub fn writeHeader(self: *Self, s: usize, d: usize, n: usize) !void {
            // Prepare a buffer for the bytes
            var buf: [@sizeOf(usize)]u8 = undefined;
            // Write the number of time steps
            try self.writeUsize(&buf, s);
            // Write the frame step
            try self.writeUsize(&buf, d);
            // Write the size of the grid
            try self.writeUsize(&buf, n);
        }
        /// Write the grid from a slice
        pub fn writeGrid(self: *Self, slice: *const Cells(F).Slice) !void {
            // Write the fields of interest as raw data
            inline for (&[_]Cells(F).Field{ .u, .v, .ro, .e }) |field| {
                const bytes = std.mem.sliceAsBytes(slice.items(field));
                _ = try self.writer.write(bytes);
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
