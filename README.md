### Notices

#### Mirrors

Repository:
- [Codeberg](https://codeberg.org/paveloom-university/Computational-Gas-Dynamics-S11-2022)
- [GitHub](https://github.com/paveloom-university/Computational-Gas-Dynamics-S11-2022)
- [GitLab](https://gitlab.com/paveloom-g/university/s11-2022/computational-gas-dynamics)

#### Prerequisites

Make sure you have installed:

- [Zig](https://ziglang.org)
- [Zigmod](https://github.com/nektro/zigmod)
- [Tracy Profiler](https://github.com/wolfpld/tracy) (optional)

#### Build

First, fetch the dependencies with `zigmod fetch`.

To build and install the library, run `zig build install`.

See `zig build --help` for more options.

#### Example

To run the example which simulates the convection process, run one of the following:

```bash
# Run with the default allocator
zig build run

# Run with the Tracy integration
zig build run -Dtracy
zig build run -Dtracy -Dtracy-depth=10

# Run with the testing allocator
zig build test
```

See `zig build --help` for more options.
