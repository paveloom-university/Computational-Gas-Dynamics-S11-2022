const std = @import("std");

const clap = @import("clap");
const lpm = @import("lpm");
const tracy = @import("tracy");

/// Type of a floating-point number
/// to use across the program
const F = f64;

// Command-line parameters
const Params = clap.parseParamsComptime(
    \\    --help
    \\      Display this help and exit.
    \\-n, --size <usize>
    \\      Size of the grid
    \\
    \\      [default: 1000].
    \\    --tau  <f64>
    \\      Time step.
    \\
    \\      [default: 1e-2].
    \\-h, --step <f64>
    \\      Grid step.
    \\
    \\      [default: 1e-2].
    \\
);

// Parse the command-line arguments
fn parseArgs() !clap.Result(clap.Help, &Params, clap.parsers.default) {
    // Parse the command-line arguments
    var diag = clap.Diagnostic{};
    var res = clap.parse(clap.Help, &Params, clap.parsers.default, .{
        .diagnostic = &diag,
    }) catch |err| {
        // Report useful error and exit
        try diag.report(std.io.getStdErr().writer(), err);
        return err;
    };
    return res;
}

// Model the process of convection
fn run(allocator: std.mem.Allocator) !void {
    // Parse the arguments
    const res = try parseArgs();
    defer res.deinit();
    // Show help if requested
    if (res.args.help) {
        const writer = std.io.getStdErr().writer();
        try writer.print("{s}\n{s}\n{s}\n\n{s}", .{
            "\u{001b}[1;32mlpm\u{001b}[m 0.1.0",
            "Pavel Sobolev <paveloom@riseup.net>",
            "Modelling convection with the large-particle method",
            "\u{001b}[0;33mPARAMETERS:\u{001b}[m\n",
        });
        return clap.help(writer, clap.Help, &Params, .{});
    }
    // Unpack the arguments
    const n = res.args.size orelse 1000;
    const tau = res.args.tau orelse 1e-2;
    const h = res.args.step orelse 1e-2;
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
        .tau = tau,
        .h = h,
        .grid = lpm.Grid(F){
            .n = n,
            .cells = cells,
        },
    };
    // Compute the evolution of the system for 1000 time steps
    model.compute(1000);
}

// Run the model with a default allocator
pub fn main() !void {
    // Mark the main function as a single frame
    tracy.frameMarkNamed("Main");
    // Prepare an allocator
    var tracy_allocator = tracy.TracyAllocator(null, 5).init(std.heap.page_allocator);
    // Run the model
    try run(tracy_allocator.allocator());
}

// Run the model with a testing allocator
test "Convection" {
    try run(std.testing.allocator);
}
