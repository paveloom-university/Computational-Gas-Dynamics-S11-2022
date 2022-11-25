const std = @import("std");

// Library package
const package = std.build.Pkg{
    .name = "lpm",
    .source = .{ .path = "src/lpm.zig" },
};

/// Build the library and the test executable
pub fn build(b: *std.build.Builder) void {
    // Define standard target options
    const target = b.standardTargetOptions(.{});
    // Define standard release options
    const mode = b.standardReleaseOptions();
    // Add the source code of the library
    const lib = b.addStaticLibrary(package.name, package.source.path);
    lib.setBuildMode(mode);
    lib.install();
    // Add the source code of the test executable
    const exe = b.addExecutable("convection", "examples/convection.zig");
    exe.setTarget(target);
    exe.setBuildMode(mode);
    exe.addPackage(package);
    exe.install();
    // Add an option to run the test executable as an integration test
    const test_step = b.step("test", "Test the convection model");
    {
        const tests = b.addTest("examples/convection.zig");
        tests.setBuildMode(mode);
        tests.addPackage(package);
        test_step.dependOn(&tests.step);
    }
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
}
