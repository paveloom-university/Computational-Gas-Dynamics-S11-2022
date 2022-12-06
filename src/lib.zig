//! An implementation of the large-particle method for the numerical
//! solution of hydrodynamic equations. Specifically, focusing on
//! solving a convection equation
//!
//! References:
//! - O.M. Belotserkovskii, Yu.M. Davydov, "The method of large particles in gas dynamics", Moscow (1982) (In Russian)

pub const Path = []const u8;

const grid = @import("grid.zig");

const Writer = @import("writer.zig").Writer;
pub const Cell = grid.Cell;
pub const Cells = grid.Cells;
pub const Grid = grid.Grid;
pub const Model = @import("model.zig").Model;
