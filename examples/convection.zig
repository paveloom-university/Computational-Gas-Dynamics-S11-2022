const std = @import("std");

const lpm = @import("lpm");
const tracy = @import("tracy");

// Model the process of convection
fn run(allocator: std.mem.Allocator) !void {
    // Initialize the grid
    var grid = try lpm.Grid(f64).init(.{
        .allocator = allocator,
        .tau = 1e-2,
        .h = 1e-2,
        .n = 1000,
    });
    defer grid.deinit();
}

// Run the model
pub fn main() !void {
    // Mark the main function as a single frame
    tracy.frameMarkNamed("Main");
    // Prepare an allocator
    var allocator = tracy.TracyAllocator(null, 5).init(std.heap.page_allocator);
    // Initialize an arena
    var arena = std.heap.ArenaAllocator.init(allocator.allocator());
    defer arena.deinit();
    // Run the model
    try run(arena.allocator());
}

// Run the model with a testing allocator
test "Convection" {
    try run(std.testing.allocator);
}
