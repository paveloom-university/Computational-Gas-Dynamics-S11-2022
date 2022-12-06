# This script animates the evolution of the system

"Check if the value of the option is the last argument"
function check_last(i)
    if i == length(ARGS)
        println("The last argument is reserved.")
        exit(1)
    end
end

"Parse the string, taking more arguments if it's quoted"
function parse_string(i)::String
    # Start from the first argument after the flag
    j = i
    # If the arguments starts with an apostrophe,
    s = if startswith(ARGS[i], "'")
        # Find the index of the argument
        # which ends with an apostrophe
        while !endswith(ARGS[j], "'")
            j += 1
        end
        # Join the arguments in one string
        # and remove the apostrophes
        chop(join(ARGS[i:j], ' '), head=1, tail=1)
    else
        # Return the next argument
        ARGS[i]
    end
    # Check for the last argument
    check_last(j)
    return s
end

# Define default values for optional arguments
POSTFIX = ""

# Parse the options
for i in eachindex(ARGS)
    # A postfix for the names of output files
    if ARGS[i] == "--postfix"
        try
            global POSTFIX = " ($(parse_string(i+1)))"
        catch
            println("Couldn't parse the value of the `--postfix` argument.")
            exit(1)
        end
    end
end

# Prepare color codes
RESET = "\e[0m"
GREEN = "\e[32m"
YELLOW = "\e[33m"

# Check for required arguments
if length(ARGS) < 1 || "--help" in ARGS
    println("""
        $(YELLOW)USAGE:$(RESET)
        { julia --project=. | ./julia.bash } scripts/coords.jl [--postfix <POSTFIX>] <DATA_PATH>

        $(YELLOW)ARGS:$(RESET)
            $(GREEN)<DATA_PATH>$(RESET)    Path to the data file

        $(YELLOW)OPTIONS:$(RESET)
            $(GREEN)--postfix <POSTFIX>$(RESET)    A postfix for the names of output files"""
    )
    exit(1)
end

# Define the path to the data file
DATA_PATH = ARGS[end]

"Padding in the output"
pad = " "^4

"Floating point type used across the script"
F = Float64

"Integer type used across the script"
I = UInt64

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots
using Plots.PlotMeasures

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(
    fontfamily="Computer Modern",
    dpi=300,
    legend=:topright,
    right_margin=20px,
    size=(420, 400)
)

# Define the paths to the data files
println(pad, "> Loading the data...")

"Data structure of the input file"
struct Data
    s::I
    n::I
    u::Vector{Matrix{F}}
    v::Vector{Matrix{F}}
    ρ::Vector{Matrix{F}}
    e::Vector{Matrix{F}}
end

# Get the fields with the vector type
vec_fields = fieldnames(Data)[3:end]

"Read the binary file"
function read_binary(path::AbstractString)::Data
    open(path, "r") do io
        # Read the number of time steps
        s = read(io, I)
        # Read the size of the matrix
        n = read(io, I)
        # Initialize the data struct
        data = Data(s, n, ntuple(_ -> [], 4)...)
        # For each time step
        for _ in 0:s
            # For each matrix field
            for field in vec_fields
                # Get the matrix
                matrix = permutedims(
                    reshape(
                        Vector{F}(reinterpret(F, read(io, sizeof(F) * n * n))),
                        (n, n)
                    )
                )
                # Push it to the appropriate vector
                push!(getfield(data, field), matrix)
            end
        end
        data
    end
end

# Read the data
data = read_binary(DATA_PATH)

println(pad, "> Plotting the initial states...")

"Plot a matrix from the specific field by index"
function plot(field::Symbol, index)
    pos = field == :ρ ? (-1.9, -2.15) : (-1.75, -1.75)
    heatmap(
        getfield(data, field)[index+1];
        annotations=(
            first(pos),
            last(pos),
            text(latexstring("$(field)_{$(index)}"), 10)
        )
    )
end

# Plot the heatmaps of the first 3 states
for field in vec_fields
    for i in 0:2
        plot(field, i)
        subscript_index = join(map(d -> '₀' + d, reverse(digits(i))))
        savefig(joinpath(PLOTS_DIR, "$(field)$(subscript_index).png"))
    end
end

# Create animations of the evolution
for field in vec_fields
    println(pad, "> Animating the evolution of `$(field)`...")
    anim = @animate for i in 0:data.s
        plot(field, i)
    end
    gif(anim, joinpath(PLOTS_DIR, "$(field).mp4"), fps=15, show_msg=false)
end

# Mark data for garbage collection
data = nothing

println()
