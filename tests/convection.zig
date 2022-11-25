const std = @import("std");

const lpm = @import("lpm");

// This integration test models the process of convection
test "Convection" {
    // Initialize the grid
    var grid = try lpm.Grid(f64).init(.{
        .allocator = std.testing.allocator,
        .tau = 1e-2,
        .h = 1e-2,
        .n = 1000,
    });
    defer grid.deinit();
}
