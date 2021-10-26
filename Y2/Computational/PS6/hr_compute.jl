using Parameters, Plots, Random, Distributions, Interpolations, Optim, LaTeXTabulars, LaTeXStrings #import the libraries we want
include("hr_model.jl") #import the functions that solve our growth model

cALPHA = 1.0
cF = 10.0
P, R  = Initialize(cALPHA,cF) #initialize primitive and results structs
@time Solve_model(P, R) #solve the model!
price_s = R.Pstar
inc_mass_s = sum(R.Mustar .- R.Mstar*P.nu)
ent_mass_s = R.Mstar
ext_mass_s = R.Mstar
dist_inc_s = R.Mustar .- R.Mstar*P.nu
dist_ent_s = R.Mstar*P.nu
# calculate labor
lab_inc_s = zeros(P.ns)
lab_ent_s = zeros(P.ns)
for i_s = 1:P.ns
    lab_inc_s[i_s] =  R.pf_n[i_s]*dist_inc_s[i_s]
    lab_ent_s[i_s] = R.pf_n[i_s]*dist_ent_s[i_s]
end
lab_inc_s = sum(lab_inc_s)
lab_ent_s = sum(lab_ent_s)
lab_frac_s = lab_ent_s/(lab_inc_s+lab_ent_s)
#dr_s = [1 0 0 0 0]

@time Solve_model_partII(P,R)
price_tv1 = R.PstarII
inc_mass_tv1 = sum(R.MustarII .- R.MstarII*P.nu)
ent_mass_tv1 = R.MstarII
ext_mass_tv1 = R.MstarII
dist_inc_tv1 = R.MustarII .- R.MstarII*P.nu
dist_ent_tv1 = R.MstarII*P.nu
# calculate labor
lab_inc_tv1 = zeros(P.ns)
lab_ent_tv1 = zeros(P.ns)
for i_s = 1:P.ns
    lab_inc_tv1[i_s] =  R.pf_nII[i_s]*dist_inc_tv1[i_s]
    lab_ent_tv1[i_s] = R.pf_nII[i_s]*dist_ent_tv1[i_s]
end
lab_inc_tv1 = sum(lab_inc_tv1)
lab_ent_tv1 = sum(lab_ent_tv1)
lab_frac_tv1 = lab_ent_tv1/(lab_inc_tv1+lab_ent_tv1)
dr_tv1 = 1*R.PrExitII


R.cALPHA = 2.0
@time Solve_model_partII(P,R)
price_tv2 = R.PstarII
inc_mass_tv2 = sum(R.MustarII .- R.MstarII*P.nu)
ent_mass_tv2 = R.MstarII
ext_mass_tv2 = R.MstarII
dist_inc_tv2 = R.MustarII .- R.MstarII*P.nu
dist_ent_tv2 = R.MstarII*P.nu
# calculate labor
lab_inc_tv2 = zeros(P.ns)
lab_ent_tv2 = zeros(P.ns)
for i_s = 1:P.ns
    lab_inc_tv2[i_s] =  R.pf_nII[i_s]*dist_inc_tv2[i_s]
    lab_ent_tv2[i_s] = R.pf_nII[i_s]*dist_ent_tv2[i_s]
end
lab_inc_tv2 = sum(lab_inc_tv2)
lab_ent_tv2 = sum(lab_ent_tv2)
lab_frac_tv2 = lab_ent_tv2/(lab_inc_tv2+lab_ent_tv2)
dr_tv2 = 1*R.PrExitII #figure out what is going on here


plot(R.s_g, R.dr_s, title = "Exit probabilities", label = "Standard")
plot!(R.s_g, dr_tv1, label = "α = 1")
plot!(R.s_g, dr_tv2, label = "α = 2")
savefig("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/Computational/PS6/exit_prob_1.png")


R.cF = 15
R.cALPHA = 1.0
@time Solve_model_partII(P,R)
@time Solve_model(P, R) #solve the model!
dr_tv1 = 1*R.PrExitII
R.cALPHA = 2.0
@time Solve_model_partII(P,R)
dr_tv2 = 1*R.PrExitII
R.dr_s[2] = 1
plot(R.s_g, R.dr_s, title = "Exit probabilities", label = "Standard")
plot!(R.s_g, dr_tv1, label = "α = 1")
plot!(R.s_g, dr_tv2, label = "α = 2")
savefig("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/Computational/PS6/exit_prob_2.png")


latex_tabular("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/Computational/PS6/tab1.tex",
    Tabular("rccc"),
    [["Variable", "Standard", "alpha=1", "alpha=2"],
    Rule(),
    ["Price Level", round(price_s; digits=2), round(price_tv1; digits=2), round(price_tv2; digits=2)],
    ["Mass of incumbents", round(inc_mass_s; digits=2), round(inc_mass_tv1; digits=2), round(inc_mass_tv2; digits=2)],
    ["Mass of exits", round(ext_mass_s; digits=2), round(ext_mass_tv1; digits=2), round(ext_mass_tv2; digits=2)],
    ["Aggregate Labor", round(lab_ent_s+lab_inc_s; digits=2), round(lab_ent_tv1+lab_inc_tv1; digits=2), round(lab_ent_tv2+lab_inc_tv2; digits=2)],
    ["Labor of Incumbents", round(lab_inc_s; digits=2), round(lab_inc_tv1; digits=2), round(lab_inc_tv2; digits=2)],
    ["Labor of Entrants", round(lab_ent_s; digits=2), round(lab_ent_tv1; digits=2), round(lab_ent_tv2; digits=2)],
    ["Frac of Labor in Entrants", round(lab_frac_s; digits=2), round(lab_frac_tv1; digits=2), round(lab_frac_tv2; digits=2)],
    Rule()])

# p=plot(G.k_grid,R.pf_v[:,1,20,1], title = string("Value(k)|K = ", G.K_grid[20]), label = "Employed, high z")
# plot!(G.k_grid,R.pf_v[:,1,20,2], label = "Unemployed, high z")
# plot!(G.k_grid,R.pf_v[:,2,20,1], label = "Employed, low z")
# plot!(G.k_grid,R.pf_v[:,2,20,2], label = "Unemployed, low z")
# savefig("./v.png")
# p=plot(G.k_grid,R.pf_k[:,1,20,1]-G.k_grid, title = string("k_prime(k) - k|K = ", G.K_grid[20]), label = "Employed, high z")
# plot!(G.k_grid,R.pf_k[:,1,20,2]-G.k_grid, label = "Unemployed, high z")
# plot!(G.k_grid,R.pf_k[:,2,20,1]-G.k_grid, label = "Employed, low z")
# plot!(G.k_grid,R.pf_k[:,2,20,2]-G.k_grid, label = "Unemployed, low z")
# savefig("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/Computational/PS5/k.png")
################################
