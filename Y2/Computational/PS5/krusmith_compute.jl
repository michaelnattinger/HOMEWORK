using Parameters, Plots, Random, Distributions, Interpolations, Optim #import the libraries we want
include("krusmith_model.jl") #import the functions that solve our growth model

P, G, S, R, idio_state, agg_state = Initialize() #initialize primitive and results structs
@time Solve_model(P, G, S, R, idio_state, agg_state) #solve the model!

p=plot(G.k_grid,R.pf_v[:,1,20,1], title = string("Value(k)|K = ", G.K_grid[20]), label = "Employed, high z")
plot!(G.k_grid,R.pf_v[:,1,20,2], label = "Unemployed, high z")
plot!(G.k_grid,R.pf_v[:,2,20,1], label = "Employed, low z")
plot!(G.k_grid,R.pf_v[:,2,20,2], label = "Unemployed, low z")
savefig("./v.png")
p=plot(G.k_grid,R.pf_k[:,1,20,1]-G.k_grid, title = string("k_prime(k) - k|K = ", G.K_grid[20]), label = "Employed, high z")
plot!(G.k_grid,R.pf_k[:,1,20,2]-G.k_grid, label = "Unemployed, high z")
plot!(G.k_grid,R.pf_k[:,2,20,1]-G.k_grid, label = "Employed, low z")
plot!(G.k_grid,R.pf_k[:,2,20,2]-G.k_grid, label = "Unemployed, low z")
savefig("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/Computational/PS5/k.png")
################################
