using Parameters, Plots #import the libraries we want
include("02Growth_model_sto.jl") #import the functions that solve our growth model

prim, res = Initialize() #initialize primitive and results structs
@time Solve_model(prim, res) #solve the model!
@unpack val_func_g, val_func_b, pol_func_g, pol_func_b = res
@unpack k_grid = prim

println("All done!")
################################
