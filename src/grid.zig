const std = @import("std");

/// A cell in the Euler grid
pub fn Cell(
    /// Type of a floating-point number
    comptime F: type,
) type {
    // Check if the provided type isn't a floating-point number type
    if (@typeInfo(F) != .Float) {
        @compileError("This type isn't allowed for floats");
    }
    return struct {
        const Self = @This();
        /// Velocity component along the X axis
        u: F,
        /// Velocity component along the X axis
        v: F,
        /// Density
        ro: F,
        /// Specific energy (energy per unit mass)
        e: F,
        /// Pressure (as a function of density and specific inner energy)
        p: F,
        /// Auxiliary velocity component along the X axis
        u_aux: F = 0,
        /// Auxiliary velocity component along the X axis
        v_aux: F = 0,
        /// Auxiliary specific energy (energy per unit mass)
        e_aux: F = 0,
        /// Create a data record from a cell (as bytes)
        pub fn record(self: *const Self) [32]u8 {
            return std.mem.toBytes([_]f64{ self.u, self.v, self.ro, self.e });
        }
    };
}

/// Cells in the Euler grid
pub fn Cells(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return std.MultiArrayList(Cell(F));
}

/// Euler grid
///
/// The returned struct emulates an N x N matrix of cells
pub fn Grid(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return struct {
        /// Size of the grid
        n: usize,
        /// The underlying data
        cells: Cells(F),
    };
}
