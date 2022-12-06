### Notices

#### Mirrors

Repository:
- [Codeberg](https://codeberg.org/paveloom-university/Computational-Gas-Dynamics-S11-2022)
- [GitHub](https://github.com/paveloom-university/Computational-Gas-Dynamics-S11-2022)
- [GitLab](https://gitlab.com/paveloom-g/university/s11-2022/computational-gas-dynamics)

#### Prerequisites

Make sure you have installed:

- [Julia](https://julialang.org)
- [Zig](https://ziglang.org)
- [Zigmod](https://github.com/nektro/zigmod)

Optionally, you might also want to install [Tracy Profiler](https://github.com/wolfpld/tracy).

#### Build

First, fetch the dependencies with `zigmod fetch`.

To build and install the library, run `zig build install`.

To run unit tests, run `zig build test`.

See `zig build --help` for more build options.

#### Example

To run the example which simulates the convection process, run one of the following:

```bash
# Run with the default allocator
zig build run -- -o res.bin
zig build run -Drelease-fast -- -o res.bin

# Run with the Tracy integration
zig build run -Dtracy -Drelease-fast -- -o res.bin
zig build run -Dtracy -Dtracy-depth=10 -Drelease-fast -- -o res.bin

# Run with the testing allocator
zig build test-run

# Run the benchmark
zig build bench
zig build bench -Drelease-fast
```

See `zig build --help` for more build options.

See `zig build run -- --help` for CLI arguments.

To animate the results, first, instantiate the project with

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Then, run one of the following

```bash
# Run without a daemon
julia --project=. scripts/animate.jl res.bin

# Run with a daemon
./julia.bash scripts/animate.jl res.bin
```
