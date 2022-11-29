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
    \\-n <usize>
    \\      Size of the grid
    \\
    \\      [default: 1000].
    \\-t, --tau <f64>
    \\      Time step.
    \\
    \\      [default: 1e-2].
    \\-h <f64>
    \\      Grid step.
    \\
    \\      [default: 1e-2].
    \\-s <usize>
    \\      Number of time steps to compute.
    \\
    \\      [default: 1000].
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

// Compute the pressure (equation of state)
inline fn pressure(
    /// Velocity component along the X axis
    u: F,
    /// Velocity component along the X axis
    v: F,
    /// Density
    ro: F,
    /// Specific energy
    e: F,
) F {
    // Compute the internal specific energy
    const eps = e - (u * u + v * v) / 2;
    // Let's assume that the gas is ideal and
    // monatomic, so the adiabatic index is
    const gamma = 3.0 / 2.0;
    // And the equation of state is
    return (gamma - 1.0) * ro * eps;
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
    const n = res.args.n orelse 1000;
    const tau = res.args.tau orelse 1e-2;
    const h = res.args.h orelse 1e-2;
    const s = res.args.s orelse 1000;
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
            // Let's assume we have neon as the gas of choice (it's monatomic)

            // The velocity is zero (the gas is at rest) [m/s]
            const u = 0;
            const v = 0;
            // The density is equal everywhere [kg/m^3]
            const ro = 0.9002;
            // Let's make the temperature change evenly (per row)
            // between -25°C and +25°C (converting to Kelvins)
            const t = -25 + @intToFloat(F, i / n) * 50 / (@intToFloat(F, n) - 1) + 273.15;
            // Specific gas constant [J/kg/K]
            // (gas constant [J/K/mol] divided by the molar mass [kg/mol])
            const r = 8.314_462_618_153_24 / 0.020_179_76;
            // The specific energy [J/kg = m^2/s^2] is
            const e = 3 / 2 * r * t;
            // The pressure [Pa = kg/m/s^2]
            const p = pressure(u, v, ro, e);
            // Initialize the cell
            cells.appendAssumeCapacity(lpm.Cell(F){
                .u = u,
                .v = v,
                .ro = ro,
                .e = e,
                .p = p,
            });
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
    model.compute(s);
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

// Benchmark the model
test "Benchmark" {
    // Set the number of iterations
    const n = 20;
    // Prepare a vector for storing the elapsed time
    var elapsed = @splat(n, @as(F, 0));
    // Initialize a timer
    var timer = try std.time.Timer.start();
    // For each iteration
    var i: usize = 0;
    while (i < n) : (i += 1) {
        // Reset the timer
        timer.reset();
        // Run the target function
        try run(std.heap.page_allocator);
        // Read the timer value
        elapsed[i] = @intToFloat(F, timer.read());
    }
    // Compute the statistics
    const min = @reduce(.Min, elapsed);
    const max = @reduce(.Max, elapsed);
    const mean = @reduce(.Add, elapsed) / n;
    const diffs = elapsed - @splat(n, mean);
    const std_dev = @sqrt(@reduce(.Add, diffs * diffs) / n);
    // Print the statistics
    std.debug.print(
        "\nIterations: {}\nMax: {:>28.15} s.\nMin: {:>28.15} s.\nMean: {:>27.15} s.\nStd. dev.: {:>22.15} s.\n",
        .{ n, min * 1e-9, max * 1e-9, mean * 1e-9, std_dev * 1e-9 },
    );
}
