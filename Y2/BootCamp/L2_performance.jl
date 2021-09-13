## show off compilation
function f(a, b)
    y = (a + 8b)^2
    return 7y
end

f(1, 2)
@code_native f(1, 2)  #view machine code generate by Julia

## Scope of variables
for i = 1:100
    println(i)
end
i # i is local within each loop, not accessible outside the loop

# globals vs not
a = 3.0 #bad global variable
@elapsed for i = 1:1000000
    global a
    a+=i
end

function timetest()
    a = 3.0 #good function-contained variable
    for i = 1:1000000
        a+=i
    end
    return a
end
@elapsed timetest() #much faster!

## showing off profiler
using Profile

function test(x::Float64)
    for i = 1:10000000
        x+=i #do a bunch of things to x
        x-=i
        x*=i
        x/=i
        x = abs(x)
        x = log(x)
    end
    x
end

Profile.clear() #empty the profile
@profile test(2.0) #profile a call of test function
Profile.print() #print profile

## Time a function
@time test(2.0)
## In-place vs out-of-place functions
# In-place functions
# modify their inputs
# Julia convention---an in-place function's name ends with `!`
# used to save memory

function f(x)
    return 2 * x
end

function f!(out, x) # example of a wrong function
    return out = 2 * x
end

y=ones(2); out = [1000;1000];
out1 = f(y)
f!(out,y); out
# The result is supposed to be 2 instead of the initial value of `out`. The f!() function works as follows:
# create a temporary variable = 2*x
# assign a local name `out` to it
# the temporary variable is not pushed to the input `out` so nothing is modified

# To correct for this error,
# Method 1: make the temporary variable 2*x global, not recommended
function f1!(out, x)
    global out = 2 * x
end
out = [1000;1000];
f1!(out,y); out
# Method 2: copy the value instead of assign the name
function f2!(out, x)
    return out[:] = 2 * x
end
out = [1000;1000];
f2!(out,y); out
# The best way is to avoid writing an in-place function and use it with caution.

## Nothing, NaN, and Missing
### Sometimes the results are absent. You can use a value "nothing" of Type "Nothing".
function f(x)
    if x >=0
        return sqrt(x)
    else
        return nothing
    end
end
println(f(-1))

# To see whether the output is nothing
isnothing(f(-1))
f(-1) == nothing

### More often, we use NaN
### NaN (not a number) is trikier. For example, Inf/Inf = NaN
exp(1000)
exp(1000)/exp(1000)
# keep in mind to use isnan() to check whether an object is NaN
isnan(NaN)
NaN == NaN

### The value missing of type Missing is used to represent missing value in a statistical sense (similar to "." in Stata).
x = [3.0, missing, 5.0, missing, missing]

using Statistics
@show mean(x)
@show mean(skipmissing(x))
@show coalesce.(x, 0.0);  # replace missing with 0.0;

## Optimal Savings
using Parameters, Plots

# Household preference: u(c) = log(c)
# Production technology: y = k^θ
# Capital depreciation rate: δ
# Budget constraint: c + k' = k^θ + (1-δ)*k

# Dynamic programming problem:
# V(k) = max_{k'} u(c) + β V(k') s.t. c + k' = k^θ + (1-δ)*k
# Solve the problem and plot the policy function k' = k'(k)
# Guess V(k'), get V(k) (function Bellman() below)
# Iterate until the function converge (function Solve_model() below)

# struct to hold model primitives (global constants)
@with_kw struct Primitives
    β::Float64 = 0.99 #discount factor
    θ::Float64 = 0.36 #production
    δ::Float64 = 0.025 #depreciation
    k_grid::Array{Float64,1} = collect(range(0.1, length = 1800, stop= 45.0)) #capital grid
    nk::Int64 = length(k_grid) #number of capital grid states
end

# struct to hold model outputs
mutable struct Results
    val_func::Array{Float64,1} #value function
    pol_func::Array{Float64,1} #policy function
end

#Bellman operator
function Bellman(prim::Primitives, res::Results)
    @unpack β, δ, θ, nk, k_grid = prim #unpack primitive structure
    v_next = zeros(nk) #preallocate next guess of value function

    for i_k = 1:nk #loop over state space
        max_util = -1e10
        k = k_grid[i_k] #value of capital
        budget = k^θ + (1-δ)*k #budget

        for i_kp = 1:nk #loop over choice set, find max {log(c) + β * res.val_func[i_kp]}
            kp = k_grid[i_kp] #value of k'
            c = budget - kp #consumption
            if c>0 #check if postiive
                val = log(c) + β * res.val_func[i_kp] #compute utility
                if val > max_util #wow new max!
                    max_util = val #reset max value
                    res.pol_func[i_k] = kp #update policy function
                end
            end
        end
        v_next[i_k] = max_util #update value function
    end
    v_next
end

#function to solve the model
function Solve_model()
    #initialize primitives and results
    prim = Primitives()
    val_func, pol_func = zeros(prim.nk), zeros(prim.nk)
    res = Results(val_func, pol_func)

    error, n = 100, 0
    while error>0.0001 #loop until convergence
        n+=1
        v_next = Bellman(prim, res) #next guess of value function
        error = maximum(abs.(v_next .- res.val_func)) #check for convergence
        res.val_func = v_next #update
        println("Current error: ", error)
    end
    println("Value function converged in ", n, " iterations")
    prim, res
end

@elapsed prim, res = Solve_model() #solve the model.
Plots.plot(prim.k_grid, res.val_func) #plot value function
