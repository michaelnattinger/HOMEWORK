## Introduction
    # Reference: [QuantEcon](https://julia.quantecon.org/)
    # [Florian Oswald's notes](https://scpo-compecon.github.io/CoursePack/)

## Execute (part of the code) in Juno
    # `Shift + Enter` will evaluate a highlighted selection or line
    # The run symbol in the menu bar (or `Ctrl+Shift+Enter`) will run the whole file

    # Example: print string "hello world!"
    println("hello world!")


## Package Environments
    # Julia’s package manager lets you set up packages used in your current code by specifying the "dependencies" (i.e., required packages).

    # An environment is a set of packages specified by a `Project.toml` (and optionally, a `Manifest.toml`). Do not delete Project.toml or Manifest.toml in your working directory!

### Download packages
    # When: first time you use it or install update

    # Example: to download package "Plots.jl", first import the `Pkg` package to your code. Then use command `Pkg.add("Plots")` to download package `"Plots.jl"`.
    import Pkg; Pkg.add("Plots")

### Use packages in your own code
using Plots

## Working with multiple scripts
# If your code is very long, you can move part of the code to a separate `.jl` file and insert the file to the main code.
include("L1_subcode.jl")

## Inserting unicode (e.g. Greek letters)
    # Julia supports the use of [unicode characters](https://docs.julialang.org/en/v1/manual/unicode-input/) such as `α` and `β` in your code

    # Unicode characters can be typed quickly using the `tab` key

    # Try creating a new code cell and typing `\alpha`, then hitting the `tab` key on your keyboard
    α=0.9;

## Data types
    # [Reference](https://docs.julialang.org/en/v1/manual/arrays/)
    # In Julia, arrays and tuples are the most important data type for working with numerical data.

### Array is mutable which means you can change its data value and modify its structure, a tuple is immutable.

### In this section we give more details on
    # creating and manipulating Julia arrays
    # fundamental array processing operations
    # basic matrix algebra
    # tuples and named tuples
    # ranges
    # nothing, missing, and unions

### One dimensional array is assumed to be a column vector. Use either semicolon or comma to connect elements.
a = [1; 2; 3]
b = [1.0, 2.0, 3.0]
    # The output tells us that the arrays are of types `Vector{Int64}` and `Vector{Float64}` respectively. Here `Int64` and `Float64` are types for the elements inferred by the compiler.

    # It is important to distinct between integers and floating point numbers because there are certain operations one type cannot perform.
        # only use integers to index the location of an entry in an array.
        a[1]
        a[1.0]

    # Display the type of `a`
    typeof(a)

### A row vector is a 2-dimensional (i.e., 1 * N matrix).
c = [1 2 3]
# The output is Array{Int64,2}. "2" means this is a 2-dim object.

### Check dimensions/sizes
size(a) # return a 1-dim turple (3,)
size(a,1) # return a number, often more useful
ndims(a) # return the dimension of a
ndims(c) # return the dimension of c

### Array vs Vector vs Matrix
    # Vector{Float64} = Array{Float64,1}
    # Matrix{Int64} = Array{Int64, 2}
    # I mentions this because people use terms interchangibly, so don't get confused when reading the documentations.
    Array{Int64, 1} == Vector{Int64}
    Array{Int64, 2} == Matrix{Int64}
    Array{Int64, 3} == Matrix{Int64}

### Creating Arrays
    # Functions that Create Arrays
    x = zeros(3)
    x = zeros(3,3)
    x = ones(3)
    x = ones(3,3)
    x = fill(1.0,3,3,3) # To return an array filled with a single value, use `fill`
    y = similar(x) # create an array of the same size

    # You can create an empty array using the `Array()` constructor
    x = Array{Float64}(undef, 3) # undef means undefined numbers

    # The printed values you see here are just garbage values. So why doing this?
    # This could be useful if you want to define the default data type in an array. For example, you cannot do this
    x = Array{Float64}(3.0, 3)

### Creating Arrays from Existing Arrays
    # Suppose you want to create a temporary variable `y = x`, this will only bind an additional name `y` to the value `x` points to.
    y = x;
    y[1] = 123;
    x

    # To create a new variable), you need to be more explicit
    # use copy()
    y = copy(x);
    y[1] = 100;
    x, y

    # use y[:] = x instead of y=x
    y[:] = x;
    y[1] = 1000;
    x, y


### Array Indexing
    x[end-1]
    x[1:end]
    y=[x x]
    y[1,:]

    # select elements that satisfy certain conditions.
    index = x .== 123;
    x[index]
    index = argmax(x);
    x[index]

### Arrays can be fairly flexible
# different types, bad for Julia's speed
x = ["foo" "bar" 1; "foo" "bar" 2]

# array of arrays
# suppose you have observations on [age gender wage] in 2 cities, you can stack observations in a city to a vector
data = [[50 1 1000.0; 60 2 6000.0],[25 1 2500.1]]
# do not need to create an additional variable "city"
# very powerful for parallelization

### tuples
x = ("foo", "bar", 1)
x[3]
x[3] = 2

## ranges, fancier indexing
x = collect(1:1:100)
y = collect(range(1, length = 100, stop = 100))
@show x[2:50]
x[75:end]

# loop over an array
languages = ["julia", "python", "stata"]
for l in languages
    println(l)
end

# iterator
for i = 1:100
    println(i)
end

# while loop
val = 1.0
tol = 0.002
while val > tol
    global val
    val = val / 2
end
val

# Conditionals
x = 1
x == 1
x == 2
x == 1 && x == 2
x == 1 || x == 2

if x == 2
    println("true")
elseif x!=2
    println("false")
end

## Making functions
function foo(x; a = 2)
    y = x^2 * a
    x, y
end

x, y = foo(2)
x, y = foo(2, a = 3)

f(x) = x^2
f(4)

# broadcasting
x = collect(1:1:100)
y = log.(x)


## Example 1: White Noise
using Plots, Distributions

dist = Normal(0,1) # Define a standard Normal distribution
n = 1000 # Sample size
ϵ = rand(dist, n) # Simulate a sample of Normal distribution using rand function
plot(1:n, ϵ) # plot(x-axis, y-axis)
histogram(ϵ) # plot distribution of the simulated sample

## Example 2: OLS
β_0 = 1.0
β_1 = 2.0
β_2 = 3.0
n = 10000
x = rand(n).*10
x2 = x.^2
ϵ = rand(dist, n)
Y = β_0 .+ β_1.*x + β_2.*x2 .+ ϵ
X = hcat(ones(n), x, x2) # horizontal catenation
β_ols = inv(X' * X) * X' * Y
