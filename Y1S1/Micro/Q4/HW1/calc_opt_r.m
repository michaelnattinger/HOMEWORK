function [ropt,popt] = calc_opt_r(grid,c,b10,b20,tol,maxiter,tune,ndraws)
options = optimoptions('fmincon','Display','off');
obj = @(r) -1*expectedprofit(r,grid,c,b10,b20,tol,maxiter,tune,ndraws);
[ropt,popt] = fmincon(obj,c,[],[],[],[],c,1,[],options);
popt = popt*-1;
end

function profit = expectedprofit(r,grid,c,b10,b20,tol,maxiter,tune,ndraws)
[b1,b2] = calc_b(grid,r,b10,b20,tol,maxiter,tune); % calculate b's
 % draw types
 draws = rand(2,ndraws);
 draws(:,2) = sqrt(draws(:,2));
 bids = 0*draws;
for i=1:ndraws
    if r<draws(1,i) 
        [~,ind] = min(abs(draws(1,i) - grid));
        bids(1,i) = b1(ind);
    end
    if r<draws(2,i) 
        [~,ind] = min(abs(draws(2,i) - grid));
        bids(2,i) = b2(ind);
    end
end
highbid = max(bids);
highbid(highbid<r) = 0;
profit = mean(highbid - c*(highbid>0),2);
end