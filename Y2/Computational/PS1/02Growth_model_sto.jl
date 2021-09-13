#keyword-enabled structure to hold model primitives
@with_kw struct Primitives
    β::Float64 = 0.99 #discount rate
    δ::Float64 = 0.025 #depreciation rate
    α::Float64 = 0.36 #capital share
    k_min::Float64 = 0.01 #capital lower bound
    k_max::Float64 = 75.0 #capital upper bound
    nk::Int64 = 3000#1000 #number of capital grid points
    k_grid::Array{Float64,1} = collect(range(k_min, length = nk, stop = k_max)) #capital grid
    Zg::Float64 = 1.25
    Zb::Float64 = 0.2
    Pgg::Float64 = 0.977
    Pbg::Float64 = 0.074
    Pbb::Float64 = 0.926
    Pgb::Float64 = 0.023
end

#structure that holds model results
mutable struct Results
    val_func_g::Array{Float64, 1} #value function
    val_func_b::Array{Float64, 1} #value function
    pol_func_g::Array{Float64, 1} #policy function
    pol_func_b::Array{Float64, 1} #policy function
end

#function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primtiives
    val_func_b = zeros(prim.nk) #initial value function guess
    pol_func_b = zeros(prim.nk) #initial policy function guess
    val_func_g = zeros(prim.nk) #initial value function guess
    pol_func_g = zeros(prim.nk) #initial policy function guess
    res = Results(val_func_b,val_func_g, pol_func_b,pol_func_g) #initialize results struct
    prim, res #return deliverables
end

#Bellman Operator
function Bellman(prim::Primitives,res::Results)
    @unpack val_func_g,val_func_b = res #unpack value function
    @unpack k_grid, β, δ, α, nk, Zg, Zb, Pgg, Pbb, Pgb, Pbg = prim #unpack model primitives
    vg_next = zeros(nk) #next guess of value function to fill
    vb_next = zeros(nk)
     #for exploiting monotonicity of policy function
    for iZ = 1:2 #1 good 2 bad
        if iZ>1
            Z = Zb
            Pb = Pbb
            Pg = Pbg
        else
            Z=Zg
            Pb = Pgb
            Pg = Pgg
        end
        choice_lower = 1
    for k_index = 1:nk
        k = k_grid[k_index] #value of k
        candidate_max = -Inf #bad candidate max
        budget = Z*k^α + (1-δ)*k #budget
        for kp_index in choice_lower:nk #loop over possible selections of k', exploiting monotonicity of policy function
            c = budget - k_grid[kp_index] #consumption given k' selection
            if c>0 #check for positivity
                val = log(c) + β*(Pg*val_func_g[kp_index] + Pb*val_func_b[kp_index] ) #compute value <-- Changes here
                if val>candidate_max #check for new max value
                    candidate_max = val #update max value
                    if iZ==1
                        res.pol_func_g[k_index] = k_grid[kp_index] #update policy function
                    else
                        res.pol_func_b[k_index] = k_grid[kp_index]
                    end
                    choice_lower = kp_index #update lowest possible choice
                end
            end
        end
        if iZ==1
            vg_next[k_index] = candidate_max #update value function
        else
            vb_next[k_index] = candidate_max
        end
    end
end
    vg_next,vb_next #return next guess of value function
end

#Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter

    while err>tol #begin iteration
        vg_next,vb_next = Bellman(prim, res) #spit out new vectors Note: I changed defn of convergence to be consistent w/ matlab
        errg = abs.(maximum(vg_next.-res.val_func_g))#/abs(vg_next[prim.nk, 1]) #reset error level
        errb = abs.(maximum(vb_next.-res.val_func_b))#/abs(vb_next[prim.nk, 1])
        err = maximum([errg errb])# errg+errb
        res.val_func_g = vg_next #update value function
        res.val_func_b = vb_next
        n+=1
    end
    println("Value function converged in ", n, " iterations.")
end

#solve the model
function Solve_model(prim::Primitives, res::Results)
    V_iterate(prim, res) #in this case, all we have to do is the value function iteration!
end
##############################################################################
