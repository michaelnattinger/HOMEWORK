using Parameters, Plots, Random, Distributions, Interpolations, Optim #import the libraries we want
include("krusmith_model.jl") #import the functions that solve our growth model

P, G, S, R, idio_state, agg_state = Initialize() #initialize primitive and results structs
@time Solve_model(P, G, S, R, idio_state, agg_state) #solve the model!
@unpack val_func_g, val_func_b, pol_func_g, pol_func_b = res
@unpack k_grid = prim

println("All done!")
################################
