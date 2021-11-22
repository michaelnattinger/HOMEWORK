function [B,J] = fmin_mn(XX,YY,B0,tol,tune)
conv = 0;
while ~conv
    B = B0 - tune*inv(hessian_mn(B0,XX))*score_mn(B0,XX,YY);
    disp(num2str(max(abs(B-B0))))
    if max(abs(B-B0))<tol
        J = LL_mn(B,XX,YY);
        conv = 1;
    else
        B0 = B;
    end
end