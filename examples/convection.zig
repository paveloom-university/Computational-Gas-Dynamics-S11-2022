const std = @import("std");

const lpm = @import("lpm");

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
    // Initialize an arena with a page allocator
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    // Run the model
    try run(arena.allocator());
}

// Run the model with a testing allocator
test "Convection" {
    try run(std.testing.allocator);
}
