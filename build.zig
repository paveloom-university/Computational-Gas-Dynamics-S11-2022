const std = @import("std");

const deps = @import("deps.zig");
const clap_pkg = deps.pkgs.clap.pkg.?;
const tracy_pkg = deps.pkgs.tracy.pkg.?;

/// Build the library and the test executable
pub fn build(b: *std.build.Builder) void {
    // Define standard target options
    const target = b.standardTargetOptions(.{});
    // Define standard release options
    const mode = b.standardReleaseOptions();
    // Add the source code of the library
    const lib = b.addStaticLibrary("lpm", "src/lib.zig");
    lib.setBuildMode(mode);
    lib.install();
    // Add the source code of the test executable
    const exe = b.addExecutable("convection", "examples/convection.zig");
    exe.setTarget(target);
    exe.setBuildMode(mode);
    exe.install();
    // Add an option to run the unit tests
    const unit_tests_step = b.step("test", "Run the unit tests");
    const unit_tests = b.addTest("src/lib.zig");
    unit_tests.setBuildMode(mode);
    unit_tests_step.dependOn(&unit_tests.step);
    // Add an option to run the test executable as an integration test
    const convection_test_step = b.step("test-run", "Test the convection model");
    const convection_test = b.addTest("examples/convection.zig");
    convection_test.setBuildMode(mode);
    convection_test.setFilter("Convection");
    convection_test_step.dependOn(&convection_test.step);
    // Add an option to run the test executable as an integration test
    const benchmark_test_step = b.step("bench", "Benchmark the convection model");
    const benchmark_test = b.addTest("examples/convection.zig");
    benchmark_test.setBuildMode(mode);
    benchmark_test.setFilter("Benchmark");
    benchmark_test_step.dependOn(&benchmark_test.step);
    // Add a run step for the executable
    const run_step = b.step("run", "Run the convection model");
    {
        const run_cmd = exe.run();
        run_cmd.step.dependOn(b.getInstallStep());
        if (b.args) |args| {
            run_cmd.addArgs(args);
        }
        run_step.dependOn(&run_cmd.step);
    }

    // Create an option to enable Tracy integration
    const tracy_enabled_option = b.option(
        bool,
        "tracy",
        "Enable Tracy integration",
    ) orelse false;
    // Create an option to override Tracy's default call stack capture depth
    const tracy_depth_option = b.option(
        c_int,
        "tracy-depth",
        "Override Tracy's default call stack capture depth",
    ) orelse 0;
    // Add these options to a separate group of options
    const tracy_options = b.addOptions();
    exe.addOptions("tracy_options", tracy_options);
    tracy_options.addOption(bool, "tracy", tracy_enabled_option);
    tracy_options.addOption(c_int, "tracy_depth", tracy_depth_option);
    // Add the options as a dependency to the Tracy package
    const new_tracy_pkg = std.build.Pkg{
        .name = "tracy",
        .source = .{ .path = tracy_pkg.source.path },
        .dependencies = &[_]std.build.Pkg{
            .{ .name = "tracy_options", .source = tracy_options.getSource() },
        },
    };
    // If the Tracy integration is enabled, link the libraries
    if (tracy_enabled_option) {
        // Gotta call this snippet until there is a nicer way
        inline for (comptime std.meta.declarations(deps.package_data)) |decl| {
            const pkg = @as(deps.Package, @field(deps.package_data, decl.name));
            for (pkg.system_libs) |item| {
                exe.linkSystemLibrary(item);
            }
            inline for (pkg.c_include_dirs) |item| {
                exe.addIncludePath(@field(deps.dirs, decl.name) ++ "/" ++ item);
            }
            inline for (pkg.c_source_files) |item| {
                exe.addCSourceFile(@field(deps.dirs, decl.name) ++ "/" ++ item, pkg.c_source_flags);
            }
        }
    }
    // Define the library package
    const new_lpm_pkg = std.build.Pkg{
        .name = "lpm",
        .source = .{ .path = "src/lib.zig" },
        .dependencies = &[_]std.build.Pkg{new_tracy_pkg},
    };
    // Add the packages
    for ([_]*std.build.LibExeObjStep{ lib, exe, unit_tests }) |step| {
        step.addPackage(new_tracy_pkg);
    }
    for ([_]*std.build.LibExeObjStep{ exe, convection_test, benchmark_test }) |step| {
        step.addPackage(clap_pkg);
        step.addPackage(new_lpm_pkg);
        // Switch to the `stage1` compiler to avoid this bug:
        // https://github.com/Hejsil/zig-clap/issues/84
        step.use_stage1 = true;
    }
}
