[ny,ny0] = size(y0);
simsim = zeros(ny,T,ny0);
for i=1:ny0
    simsim(:,1,i) = y0(:,i); % Initial values
end

% solve for equations we can use to project exactly into the future
for i=1:ny0 % for each sim corresponding to each initial value...
for tt=2:T  % for each point in time...
simsim(1,tt,i) = param.pz * simsim(1,tt-1,i)^(param.palpha) + (1 - param.pdelta)*simsim(1,tt-1,i) - simsim(2,tt-1,i);
simsim(1,tt,i) = max([0.001 simsim(1,tt,i)]); % must be positive
simsim(2,tt,i) = param.pbeta*simsim(2,tt-1,i)*(1 - param.pdelta + param.palpha*param.pz * simsim(1,tt,i)^(param.palpha - 1));
simsim(2,tt,i) = max([0.001 simsim(2,tt,i)]); % must be positive
end
end