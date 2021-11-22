function quad = LL_quad(XX,YY,ZZ,params,nodes1,weights1,nodes2,weights2)
[T,nx] = size(XX);
quad = 0;
sig = 1/(1-params.rho)^2;
sig = sqrt(sig);

for i_t=1:T % for each point in time
if YY(i_t,1) 
    prob = normcdf(-params.alpha0 - XX(i_t,:)*params.beta - ZZ(i_t,1)*params.gamma);
    quad = quad + log(prob);
elseif YY(i_t,2)
    ub = params.alpha0 + XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma;
    PHI = -params.alpha1 - XX(i_t,:)*params.beta - ZZ(i_t,2)*params.gamma;
    prob = int_1d(ub,PHI,sig,params.rho,nodes1,weights1);
    quad = quad + log(prob);
elseif YY(i_t,3)
    ubl = params.alpha0 +  XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma;
    ubr = params.alpha1 +  XX(i_t,:)*params.beta + ZZ(i_t,2)*params.gamma;
    PHI = -params.alpha2 - XX(i_t,:)*params.beta - ZZ(i_t,3)*params.gamma;
    prob = int_2d(ubl,ubr,PHI,sig,params.rho,nodes2,weights2);
    quad = quad + log(prob);
else 
    ubl = params.alpha0 +  XX(i_t,:)*params.beta + ZZ(i_t,1)*params.gamma;
    ubr = params.alpha1 +  XX(i_t,:)*params.beta + ZZ(i_t,2)*params.gamma;
    PHI = params.alpha2 + XX(i_t,:)*params.beta + ZZ(i_t,3)*params.gamma;
    prob = int_2d(ubl,ubr,PHI,sig,params.rho,nodes2,weights2);
    quad = quad + log(prob);
end
%if isinf(abs(GHK)); GHK = -10e10; end
end
end

function prob = int_1d(ub,PHI,sig,rho,nodes1,weights1)
t_nodes = log(nodes1)+ub;
j_nodes = 1./nodes1;
prob = sum(normcdf(PHI - rho.*t_nodes).*(normpdf(t_nodes/sig)/sig).*weights1.*j_nodes);
end

function prob = int_2d(ubl,ubr,PHI,sig,rho,nodes2,weights2)
t_nodes = log(nodes2) + [ubl ubr];
j_nodes = 1./nodes2;
prob = sum(normcdf(PHI - rho.*t_nodes(:,2)).*normpdf(t_nodes(:,2) - rho.*t_nodes(:,1)).*(normpdf(t_nodes(:,1)./sig)./sig).*weights2.*prod(j_nodes,2));
end