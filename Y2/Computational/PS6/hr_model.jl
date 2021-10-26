#keyword-enabled structure to hold model primitives
@with_kw struct Params
    cBET::Float64 = 0.8
    cTHETA::Float64 = 0.64
    cA::Float64 = 1/200
    #cF::Float64 = 10 #going to vary eventually
    cE::Float64 = 5
    #cALPHA::Float64 = 1 #going to vary eventually
    ns::Int64 = 5
    s_grid::Array{Float64,1} = [3.98*10^(-4), 3.58, 6.82, 12.18, 18.79]
    F::Array{Float64,2} = [0.6598 0.26 0.0416  0.0331  0.0055; 0.1997 0.7201 0.0420 0.0326 0.0056; 0.2 0.2 0.5555 0.0344 0.0101; 0.2 0.2 0.2502 0.3397 0.0101; 0.2 0.2 0.25 0.34 0.01]
    nu::Array{Float64,1} = [0.37, 0.4631, 0.1102, 0.0504, 0.0063]
end

mutable struct Results
    pf_v::Array{Float64,1}
    pf_x::Array{Int64,1}
    pf_n::Array{Float64,1}
    pf_pr::Array{Float64,1}
    Pstar::Float64
    Mstar::Float64
    Mustar::Array{Float64,1}
    EC::Float64
    LMC::Float64
    pf_vII::Array{Float64,2}
    pf_UII::Array{Float64,1}
    pf_sigII::Array{Float64,1}
    pf_nII::Array{Float64,1}
    pf_prII::Array{Float64,1}
    PstarII::Float64
    MstarII::Float64
    MustarII::Array{Float64,1}
    ECII::Float64
    LMCII::Float64
    cALPHA::Float64
    cF::Float64
    PrExitII::Array{Float64,1}
    s_g::Array{Float64,1}
    dr_s::Array{Float64,1}
end


function Initialize(input1::Float64,input2::Float64)
    P = Params() #initialize primtiives
    pf_x = zeros(P.ns) #initial value function guess i_k, i_eps, i_K, i_z
    pf_v = zeros(P.ns) #initial policy function guess
    pf_n = zeros(P.ns)
    pf_pr = zeros(P.ns)
    Pstar = 0.738
    Mstar = 2.6386
    Mustar = zeros(P.ns)
    EC = 1
    LMC = 1

    pf_vII = zeros(P.ns,2)
    pf_UII = zeros(P.ns)
    pf_sigII = zeros(P.ns)
    pf_nII = zeros(P.ns)
    pf_prII = zeros(P.ns)
    PstarII = 0.6908
    MstarII = 3.355
    MustarII = zeros(P.ns)
    PrExitII = zeros(P.ns)
    ECII = 1
    LMCII = 1
    cALPHA = input1
    cF = input2
    s_g = (1:P.ns)
    dr_s = zeros(P.ns)
    dr_s[1] = 1


    R = Results(pf_v,pf_x,pf_n,pf_pr,Pstar,Mstar,Mustar,EC,LMC, pf_vII, pf_UII, pf_sigII, pf_nII,pf_prII, PstarII, MstarII, MustarII, ECII, LMCII,cALPHA,cF,PrExitII,s_g,dr_s) #initialize results struct
    P, R
end

function Bellman(P::Params,  R::Results)
    @unpack cBET, cTHETA, cA, cE, ns, s_grid, F, nu = P

    @unpack pf_v,pf_x,pf_n,pf_pr,Pstar,Mstar,Mustar,EC,LMC,cF= R

    pf_x_up = zeros(ns)
    pf_v_up = zeros(ns)
    pf_n_up = zeros(ns)
    pf_pr_up = zeros(ns)
    contin = 0.0
    for i_s = 1:ns
        contin = 0.0
        for i_sp = 1:ns
        contin = contin +  F[i_s,i_sp]*pf_v[i_sp]
        end
        if contin>=0
            pf_pr_up[i_s],pf_n_up[i_s] = profits(Pstar,s_grid[i_s],cF,cTHETA)
            pf_v_up[i_s] = pf_pr_up[i_s] + cBET * contin
            pf_x_up[i_s] = 0
        else
            pf_pr_up[i_s],pf_n_up[i_s] = profits(Pstar,s_grid[i_s],cF,cTHETA)
            pf_x_up[i_s] = 1
            pf_v_up[i_s] = pf_pr_up[i_s]
        end
    end
    return pf_x_up, pf_v_up, pf_n_up, pf_pr_up
end

function profits(Pstar::Float64,state::Float64, cF::Float64, cTHETA::Float64)
    n = (Pstar*state*cTHETA)^(1/(1-cTHETA))
    n = maximum([n 0])
    prf = Pstar*state*n^(cTHETA) - n - Pstar*cF
    return prf, n
end

#Value function iteration
function V_iterate(P::Params, R::Results; tol::Float64 = 1e-9, err::Float64 = 100.0)
    n = 0 #counter
    while err>tol #begin iteration
        pf_x_up, pf_v_up, pf_n_up, pf_pr_up = Bellman(P,R) #spit out new vectors
        err = abs.(maximum(pf_v_up.-R.pf_v))#/abs(vg_next[prim.nk, 1]) #reset error level
        R.pf_v = pf_v_up #update value function
        R.pf_x = pf_x_up
        R.pf_n = pf_n_up
        R.pf_pr = pf_pr_up
        n+=1
    end
    #println("Value function converged in ", n, " iterations.")
end

#solve the model
function Solve_model(P::Params, R::Results)
    # solve for p*
    #Pstar = 1.0
    tol = 1e-6
    err = 9.0
    tune = 0.0001 # amount to increase/decrease p by
    while err>tol
        V_iterate(P,R) # VFI to get value functions
        R.EC = 0.0
        for i_s = 1:P.ns
            R.EC = R.EC + (P.nu[i_s]*R.pf_v[i_s])/R.Pstar
        end
        R.EC = R.EC - P.cE
        err = abs(R.EC)
        #println("Pstar, ", R.Pstar)
        #println("EC, ", R.EC)
        if (R.EC>0 && err>tol )
            R.Pstar = R.Pstar *(1-tune)
        elseif (R.EC<0 && err>tol )
            R.Pstar = R.Pstar *(1+tune)
        end
    end

    Mtune = 0.001
    Merr = 9.0
    Mtol = 1e-4
    while Merr>Mtol
    muerr = 9.0
    mutol = 1e-12
    mu_init = P.nu
    mu_up = 0*P.nu
    while muerr>mutol
        muerr = 0.0
        mu_up = R.Mstar*P.nu # my notation is just slightly different than what they have but results are same
        for i_s = 2:P.ns # note: starts from 2 bc mass at first index drops out
            for i_sp = 1:P.ns
                mu_up[i_sp] = mu_up[i_sp] + P.F[i_s,i_sp]*mu_init[i_s] # mass of old firms
                #mu_up[i_sp] = mu_up[i_sp] + R.Mstar* P.F[i_s,i_sp]*P.nu[i_s] # mass of new firms
            end
        end

        muerr = maximum(abs.(mu_up .- mu_init))
        mu_init = mu_up
    end
    R.Mustar = mu_up
    LD = 0.0
    PR = 0.0
    for i_s = 1:P.ns # labor demand of firm, and total profits
        LD = LD + R.pf_n[i_s]*R.Mustar[i_s] #+ R.pf_n[i_s]*P.nu[i_s]
        PR = PR + R.pf_pr[i_s]*R.Mustar[i_s] #+ R.pf_pr[i_s]*P.nu[i_s]
    end
    LS = (1/P.cA) - PR
    R.LMC = LD - LS
    Merr = abs(R.LMC)
    if (Merr>Mtol && R.LMC>0.0)
        R.Mstar = R.Mstar * (1-Mtune)
    elseif (Merr>Mtol && R.LMC<0.0)
        R.Mstar = R.Mstar * (1+Mtune)
    else
        println("Mstar, ", R.Mstar)
        println("LMC, ", R.LMC)
    end
    end

end



function Solve_model_partII(P::Params, R::Results)
    tol = 1e-4
    err = 9.0
    tune = 0.001 # amount to increase/decrease p by
    while err>tol
        V_iterateII(P,R) # VFI to get value functions
        R.ECII = 0.0
        for i_s = 1:P.ns
            R.ECII = R.ECII + (P.nu[i_s]*R.pf_UII[i_s])/R.PstarII
        end
        R.ECII = R.ECII - P.cE
        err = abs(R.ECII)
        if (R.ECII>0 && err>tol )
            R.PstarII = R.PstarII *(1-tune)
        elseif (R.ECII<0 && err>tol )
            R.PstarII = R.PstarII *(1+tune)
        end
    end

    Mtune = 0.001
    Merr = 9.0
    Mtol = 1e-4
    while Merr>Mtol
    muerr = 9.0
    mutol = 1e-12
    mu_init = P.nu
    mu_up = 0*P.nu
    while muerr>mutol
        muerr = 0.0
        mu_up = R.MstarII*P.nu # my notation is just slightly different than what they have but results are same
        for i_s = 1:P.ns # note: starts from 2 bc mass at first index drops out
            for i_sp = 1:P.ns
                if (R.PrExitII[i_s]<1*10^(-12) && i_s==1)
                    println("error! No probability of exit in worst state!")
                end
                mu_up[i_sp] = mu_up[i_sp] + P.F[i_s,i_sp]*mu_init[i_s]*(1-R.PrExitII[i_s]) # mass of old firms
            end
        end
        muerr = maximum(abs.(mu_up .- mu_init))
        mu_init = mu_up
    end
    R.MustarII = mu_up
    LD = 0.0
    PR = 0.0
    for i_s = 1:P.ns # labor demand of firm, and total profits
        LD = LD + R.pf_nII[i_s]*R.MustarII[i_s] #+ R.pf_n[i_s]*P.nu[i_s]
        PR = PR + R.pf_prII[i_s]*R.MustarII[i_s] #+ R.pf_pr[i_s]*P.nu[i_s]
    end
    LS = (1/P.cA) - PR
    R.LMCII = LD - LS
    Merr = abs(R.LMCII)
    if (Merr>Mtol && R.LMCII>0.0)
        R.MstarII = R.MstarII * (1-Mtune)
    elseif (Merr>Mtol && R.LMC<0.0)
        R.MstarII = R.MstarII * (1+Mtune)
    else
        println("MstarII, ", R.MstarII)
        println("LMCII, ", R.LMCII)
    end
    end
end

function V_iterateII(P::Params, R::Results; tol::Float64 = 1e-9, err::Float64 = 100.0)
    n = 0 #counter
    while err>tol #begin iteration
        pf_U_up,  pf_v_up, pf_n_up, pf_pr_up = BellmanII(P,R) #spit out new vectors
        err = abs.(maximum(pf_U_up.-R.pf_UII))#/abs(vg_next[prim.nk, 1]) #reset error level
        R.pf_vII = pf_v_up #update value function
        R.pf_nII = pf_n_up
        R.pf_prII = pf_pr_up
        R.pf_UII = pf_U_up
        n+=1
    end
end


function BellmanII(P::Params,  R::Results)
    @unpack cBET, cTHETA, cA, cE, ns, s_grid, F, nu = P
    @unpack pf_UII,pf_vII,pf_nII,pf_prII,PstarII,MstarII,MustarII,ECII,LMCII, cF, cALPHA= R

    pf_v_up = zeros(ns,2)
    pf_n_up = zeros(ns)
    pf_pr_up = zeros(ns)
    pf_U_up = zeros(ns)
    contin = 0.0
    for i_s = 1:ns
        contin = 0.0
        for i_sp = 1:ns
        contin = contin +  F[i_s,i_sp]*pf_UII[i_sp]
        end

            pf_pr_up[i_s],pf_n_up[i_s] = profits(PstarII,s_grid[i_s],cF,cTHETA)
            pf_v_up[i_s,1] = pf_pr_up[i_s] + cBET * contin
            pf_v_up[i_s,2] = pf_pr_up[i_s]
            vMax = cALPHA*maximum(pf_v_up[i_s,:], dims=1)/6
            pf_U_up[i_s] = 0.5772156649/cALPHA + (1/cALPHA)*(vMax[1]+log(exp(cALPHA*pf_v_up[i_s,1] -vMax[1]) + exp(cALPHA*pf_v_up[i_s,2]-vMax[1])))
            R.PrExitII[i_s] = exp(cALPHA*pf_v_up[i_s,2])/(exp(cALPHA*pf_v_up[i_s,2])+exp(cALPHA*pf_v_up[i_s,1]))
    end
    return pf_U_up, pf_v_up, pf_n_up, pf_pr_up
end
##############################################################################
