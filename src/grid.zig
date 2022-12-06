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
        /// Auxiliary density
        ro_aux: F = 0,
        /// Auxiliary specific energy (energy per unit mass)
        e_aux: F = 0,
    };
}

/// Direction
pub const Dir = enum {
    top,
    bottom,
    left,
    right,
};

/// Cells in the Euler grid
pub fn Cells(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return std.MultiArrayList(Cell(F));
}

/// Euler grid
///
/// The returned struct emulates an N x N matrix of cells.
/// Note that the index grows left-to-right, bottom-to-top.
pub fn Grid(
    /// Type of a floating-point number
    comptime F: type,
) type {
    return struct {
        const Self = @This();
        /// Size of the grid
        n: usize,
        /// The underlying data
        cells: Cells(F),
        /// Check if the cell is near the edge of the grid
        pub inline fn nearEdge(self: *Self, comptime dir: Dir, index: usize) bool {
            return switch (dir) {
                .top => index / self.n == self.n - 1,
                .bottom => index / self.n == 0,
                .left => index % self.n == 0,
                .right => index % self.n == self.n - 1,
            };
        }
        /// Get a value from the neighbour cell
        pub inline fn neighbour(self: *Self, comptime dir: Dir, comptime normal: bool, array: []F, index: usize) F {
            return if (self.nearEdge(dir, index))
                // If requested a normal value, change the sign
                //
                // This should only apply to the velocities
                if (normal) -array[index] else array[index]
            else switch (dir) {
                .top => array[index + self.n],
                .bottom => array[index - self.n],
                .left => array[index - 1],
                .right => array[index + 1],
            };
        }
        /// Compute a value at the border of the cell
        pub inline fn border(self: *Self, comptime dir: Dir, comptime normal: bool, array: []F, index: usize) F {
            return (array[index] + self.neighbour(dir, normal, array, index)) / 2;
        }
    };
}
