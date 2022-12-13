const std = @import("std");

const clap = @import("clap");
const lpm = @import("lpm");
const tracy = @import("tracy");

// We use evented IO for multithreading
//
// Unfortunately, this makes things much slower
// for single-threaded scenarios. You might want
// to comment this line if you use one thread
// (which is the default!)
pub const io_mode = .evented;

/// Type of a floating-point number
/// to use across the program
const F = f64;

// Default values of the parameters
//
// Characteristic value of the velocity and the time steps
// here should follow the Courant–Friedrichs–Lewy condition
const n_default = 100;
const tau_default = 1e-8;
const h_default = 1e-4;
const s_default = 100000;
const d_default = 1000;
const seed_default = 0;
const phi_default = 6e2;
const j_default = 1;

/// Command-line arguments
const Args = struct {
    n: usize = n_default,
    tau: F = tau_default,
    h: F = h_default,
    s: usize = s_default,
    d: usize = d_default,
    seed: u64 = seed_default,
    phi: F = phi_default,
    j: usize = j_default,
    path: []const u8,
};

/// Command-line parameters
const Params = clap.parseParamsComptime(std.fmt.comptimePrint(
    \\    --help
    \\      Display this help and exit.
    \\
    \\-n <usize>
    \\      Size of the grid
    \\
    \\      [default: {}].
    \\
    \\-t, --tau <f64>
    \\      Time step [s].
    \\
    \\      [default: {}].
    \\
    \\-h <f64>
    \\      Grid step [m].
    \\
    \\      [default: {}].
    \\
    \\-s <usize>
    \\      Number of time steps to compute.
    \\
    \\      [default: {}].
    \\
    \\-d <usize>
    \\      Save results every `d` frames.
    \\
    \\      [default: {}].
    \\
    \\    --seed <u64>
    \\      Seed of the random number generator.
    \\
    \\      [default: {}].
    \\
    \\    --phi <f64>
    \\      Value of the gravitational potential [m^2/s^2].
    \\
    \\      [default: {}].
    \\
    \\-j <usize>
    \\      Number of threads to use.
    \\
    \\      [default: {}].
    \\
    \\-o, --output <str>
    \\      A path to the output file.
    \\
, .{
    n_default,
    tau_default,
    h_default,
    s_default,
    d_default,
    seed_default,
    phi_default,
    j_default,
}));

// Parse the command-line arguments
fn parseArgs() !Args {
    // Parse the command-line arguments
    var diag = clap.Diagnostic{};
    var res = clap.parse(clap.Help, &Params, clap.parsers.default, .{
        .diagnostic = &diag,
    }) catch |err| {
        // Report an error and exit
        try diag.report(std.io.getStdErr().writer(), err);
        return err;
    };
    defer res.deinit();
    // Show help if requested
    if (res.args.help) {
        const writer = std.io.getStdErr().writer();
        try writer.print("{s}\n{s}\n{s}\n\n{s}", .{
            "\u{001b}[1;32mlpm\u{001b}[m 0.1.0",
            "Pavel Sobolev <paveloom@riseup.net>",
            "Modelling convection with the large-particle method",
            "\u{001b}[0;33mOPTIONS:\u{001b}[m\n",
        });
        try clap.help(writer, clap.Help, &Params, .{});
        std.os.exit(0);
    }
    // Unpack the path
    const path = res.args.output orelse {
        const stderr = std.io.getStdErr().writer();
        try stderr.writeAll("A path to the output file is required.\n");
        std.os.exit(1);
    };
    // Return the arguments
    return Args{
        .n = res.args.n orelse n_default,
        .tau = res.args.tau orelse tau_default,
        .h = res.args.h orelse h_default,
        .s = res.args.s orelse s_default,
        .d = res.args.d orelse d_default,
        .seed = res.args.seed orelse seed_default,
        .phi = res.args.phi orelse phi_default,
        .j = res.args.j orelse j_default,
        .path = path,
    };
}

// Equation of state
fn eqs(
    /// Velocity component along the X axis
    u: F,
    /// Velocity component along the Y axis
    v: F,
    /// Density
    ro: F,
    /// Specific energy
    e: F,
) F {
    // Compute the internal specific energy [J/kg = m^2/s^2]
    const eps = e - (u * u + v * v) / 2;
    // Let's assume that the gas is ideal and
    // monatomic, so the adiabatic index is
    const gamma = 5.0 / 3.0;
    // And the equation of state is [Pa = kg/m/s^2]
    return (gamma - 1.0) * ro * eps;
}

// Model the process of convection
fn run(allocator: std.mem.Allocator, args: *const Args) !void {
    // Prepare cells
    var cells = cells: {
        // Initialize an empty list
        var cells = lpm.Cells(F){};
        // Compute the number of cells
        const m = args.n * args.n;
        // Prepare space for storing N x N cells
        try cells.ensureTotalCapacity(allocator, m);
        // Prepare a random number generator
        var rng = std.rand.DefaultPrng.init(args.seed).random();
        // Initialize each cell
        var i: usize = 0;
        while (i < m) : (i += 1) {
            // Let's assume we have neon as the gas of choice (it's monatomic)

            // The velocity components fluctuate around zero [m/s]
            const u = rng.floatNorm(F);
            const v = rng.floatNorm(F);
            // The density is equal everywhere [kg/m^3]
            const ro = 0.9002;
            // Let's make the temperature change evenly (per row)
            // between -25°C and -15°C (converting to Kelvins)
            const t = -15 -
                @intToFloat(F, i / args.n) *
                10 / @intToFloat(F, args.n - 1) +
                273.15;
            // Specific gas constant [J/kg/K]
            // (gas constant [J/K/mol] divided by the molar mass [kg/mol])
            const r = 8.314_462_618_153_24 / 0.020_179_76;
            // The specific energy [J/kg = m^2/s^2] is
            const e = 3.0 / 2.0 * r * t;
            // The pressure [Pa = kg/m/s^2]
            const p = eqs(u, v, ro, e);
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
    // Initialize the model
    var model = try lpm.Model(F).init(.{
        .allocator = allocator,
        .tau = args.tau,
        .h = args.h,
        .grid = lpm.Grid(F){
            .n = args.n,
            .cells = cells,
        },
        .phi = args.phi,
        .eqs = eqs,
        .path = args.path,
        .threads = args.j,
    });
    defer model.deinit();
    // Compute the evolution of the system for `s` time steps
    try model.compute(args.s, args.d);
}

// Run the model with a default allocator
pub fn main() !void {
    // Mark the main function as a single frame
    tracy.frameMarkNamed("Main");
    // Prepare an allocator
    var tracy_allocator = tracy.TracyAllocator(null, 5).init(std.heap.page_allocator);
    // Parse the arguments
    const args = try parseArgs();
    // Run the model
    try run(tracy_allocator.allocator(), &args);
}

// Run the model with a testing allocator
test "Convection" {
    // Prepare test arguments
    const args = Args{
        .path = "test.bin",
        .s = 1000,
        .d = 10,
    };
    // Run the model
    try run(std.testing.allocator, &args);
}

// Benchmark the model
test "Benchmark" {
    // Set the number of iterations
    const n = 20;
    // Prepare test arguments
    const args = Args{
        .path = "bench.bin",
        .s = 1000,
        .d = 10,
    };
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
        try run(std.heap.page_allocator, &args);
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
        "\nIterations: {}\ns:          {}\nd:          {}\nMax: {:>28.15} s.\nMin: {:>28.15} s.\nMean: {:>27.15} s.\nStd. dev.: {:>22.15} s.\n",
        .{ n, args.s, args.d, min * 1e-9, max * 1e-9, mean * 1e-9, std_dev * 1e-9 },
    );
}
