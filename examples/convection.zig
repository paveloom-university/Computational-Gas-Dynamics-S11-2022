const std = @import("std");

const lpm = @import("lpm");
const tracy = @import("tracy");

/// Type of a floating-point number
/// to use across the program
const F = f64;

// Model the process of convection
fn run(allocator: std.mem.Allocator) !void {
    // Set the size of the grid
    const n = 1000;
    // Prepare cells
    var cells = cells: {
        // Initialize an empty list
        var cells = lpm.Cells(F){};
        // Compute the number of cells
        const m = n * n;
        // Prepare space for storing N x N cells
        try cells.ensureTotalCapacity(allocator, m);
        // Initialize each cell
        var i: usize = 0;
        while (i < m) : (i += 1) {
            cells.appendAssumeCapacity(lpm.Cell(F){});
        }
        break :cells cells;
    };
    defer cells.deinit(allocator);
    // Initialize the model
    var model = lpm.Model(F){
        .tau = 1e-2,
        .h = 1e-2,
        .grid = lpm.Grid(F){
            .n = 1000,
            .cells = cells,
        },
    };
    // Compute the evolution of the system for 1000 time steps
    model.compute(1000);
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
