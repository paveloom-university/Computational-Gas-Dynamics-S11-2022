const std = @import("std");

// Library package
const package = std.build.Pkg{
    .name = "lpm",
    .source = .{ .path = "src/lpm.zig" },
};

/// Build the library
pub fn build(b: *std.build.Builder) void {
    // Define standard release options
    const mode = b.standardReleaseOptions();
    // Add the source code of the library
    const lib = b.addStaticLibrary(package.name, package.source.path);
    lib.setBuildMode(mode);
    lib.install();
    // Add options for running library tests
    const test_step = b.step("test", "Run the convection test");
    // Add the integration test for modelling convection
    {
        const tests = b.addTest("tests/convection.zig");
        tests.setBuildMode(mode);
        tests.addPackage(package);
        test_step.dependOn(&tests.step);
    }
}
