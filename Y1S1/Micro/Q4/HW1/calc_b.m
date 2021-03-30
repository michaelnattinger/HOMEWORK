function [b1,b2] = calc_b(grid,r,b10,b20,tol,maxiter,tune)
b1 = 0*b10;
b2 = 0*b20;
b1p= b1(2:end);
b2p=b1p;
iter = 1;
diff = 999;
ng = length(grid);
while (diff>tol)&&(iter<maxiter)
    for i=2:ng % Calculate derivatives
        [~,ind1] = min(b20-b10(i));
        [~,ind2] = min(b10-b20(i));
        b1p(i-1) = (grid(ind1) - b10(i))/grid(i);
        b2p(i-1) = 2*(grid(ind2) - b20(i))/grid(i);
    end
    b1 = numint(b1p,grid,r);
    b2 = numint(b2p,grid,r);
    diff = sum(abs(b2-b20)+abs(b1-b10)); 
    iter = iter+1;
    b10 = tune*b1+(1-tune)*b10;
    b20 = tune*b2+(1-tune)*b20;
end
end

function int = numint(fp,grid,r)
int = cumsum([r fp.*(grid(2:end)-grid(1:end-1))]);
end