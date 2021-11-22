using XLSX, Parameters, Plots, Random, Distributions, Interpolations, Optim #import the libraries we want
include("probit_funcs.jl") #import the functions that solve our growth model
# read-in data from excel
xlsdata = zeros(16401,18)
xf = XLSX.readxlsx("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/Computational/PS8/data.xlsx")
sh = xf["Sheet1"]
data = sh["A1:R16401"]
data = convert(Array{Float64,2},data)
D, R = Initialize(data)
y = Lambda(0.0)
LL = LogLiki(D,R)
