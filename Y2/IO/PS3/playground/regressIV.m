

function f = regressIV(y,XI,X,Z,fid);
%This function estimates an IV regression and standard errors.
%It is a bit slow but seems to work fine
%Inputs:
%Y - the dependent variable
%XI - other endogenous varibales that require instrumenting
%X - exogenous variables used in the second stage regression
%Z - instrumets used in the first stage regression
%fid - a dummy variable =1 to give output =0 otherwise
%Important Outputs:
%f.b - estimated coefficients
%f.seb - standard errors
% Written by Abe Dunn, 6/14/03


Xhat=Z*(Z'*Z)^(-1)*Z'*XI;
Xhat=[Xhat X];
X=[XI X];
f.Xhat=Xhat;
f.b=(Xhat'*Xhat)^(-1)*Xhat'*y;  %Greene pg374
f.e=y-X*f.b;
f.N=size(X,1);
f.K=size(X,2);
f.se2=(f.e'*f.e)/(f.N-f.K);
f.varb=f.se2*(Xhat'*Xhat)^(-1);
f.seb=f.varb.^(1/2);
f.R2=1-f.e'*f.e/(y'*y-f.N*mean(y)^2);
f.adjR2=1-(f.N-1)/(f.N-f.K)*(1-f.R2);

if fid
   fprintf(fid,'Estimates: \n');
   for i=1:f.K
   	fprintf(fid,'\n  b%d:  %10.8f',i,f.b(i));
      fprintf(fid,'\n   (%10.8f) \n', f.seb(i,i));
   end;
   fprintf('\n\n R-squared: %3.3f  \n adjR-squared: %3.3f \n\n',f.R2, f.adjR2);
end;

