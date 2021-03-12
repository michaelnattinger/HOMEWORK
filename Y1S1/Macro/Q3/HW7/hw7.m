clear; close all; clc
pbeta = 0.99;
delk = 1e-4;
delb = 1e-4;
Bmax = 1;
k = (delk:delk:1-delk)';
B = (delb:delb:Bmax)';
PI = [0.8 0.2; 0.2 0.8];
a = [0.5 1];
na = length(a);
nk = length(k);
nb = length(B);
Q0 = 1*ones(nk,na,nb); % initial guess of Q
tune = 0.99;
tol = 1e-3;
diff = 99;
maxiter = 1e4;
iter = 1;
Qmin = Q0;
Bch = Q0;
Kch = Q0;
while (diff>tol)&&(iter<maxiter)
    if mod(iter,100)==0; disp(['Iter ' num2str(iter) ', diff = ' num2str(diff)]); end
    Q=0*Q0;
    Kt = Q;
    for ai = 1:na
        for ki = 1:nk
            for bi = 1:nb
            Kt = (1/(Q0(ki,ai,bi)*(1-pbeta*(Qmin(ki,ai,bi)/Q0(ki,ai,bi)))))*((a(ai)+Q0(ki,ai,bi)).*k(ki) -(1/pbeta)*B(bi));
            [~,iKt] = min(abs(Kt - k));
            Bt = pbeta*Kt*Qmin(ki,ai,bi);
            [~,iBt] = min(abs(Bt - B));
            Q(ki,ai,bi) = pbeta*(PI(ai,1)*Q0(iKt,1,iBt)+PI(ai,2)*Q0(iKt,2,iBt))+ (pbeta/2)*(1-k(iKt));
            Qmin(ki,ai,bi) = min([Q0(iKt,1,iBt) Q0(iKt,2,iBt)]);
            Bch(ki,ai,bi) = B(iBt);
            Kch(ki,ai,bi) = k(iKt);
            %if Kt<0; error('negative kt!'); end
            end
        end
    end
    diff = sum(abs(Q(:)-Q0(:)));
    iter = iter + 1;
    Q0 = Q*(1-tune) + Q0*(tune);
end