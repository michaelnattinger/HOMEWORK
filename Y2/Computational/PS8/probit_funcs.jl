mutable struct Results
    BB::Array{Float64,2}
end

mutable struct Data
    YY::Array{Float64,2}
    XX::Array{Float64,2}
    nx::Int64
    T::Int64
end

function Initialize(xlsdata:: Array{Float64,2})
    Y = xlsdata[:,1:1]
    X = xlsdata[:,2:end]
    nx = size(X,2)
    T = size(X,1)
    D = Data(Y,X,nx,T)
    B0 = zeros(nx,1)
    B0[1,1] = -1.0
    R = Results(B0)
    D,R
end

function LogLiki(D::Data,R::Results)
    sum = 0.0
for tt=1:D.T
    Prob = Lambda(D.XX[tt:tt,:]*R.BB)
    sum = sum + log(Prob^D.YY[tt,1] * (1-Prob)^(1-D.YY[tt,1]))
end
    sum
end

function Lambda(x::Float64)
    y = exp(x)/(1+exp(x))
    y
end

function Lambda(x::Array{Float64,2})
    x = x[1,1]
    y = exp(x)/(1+exp(x))
    y
end
