function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g2
%

if T_flag
    T = noadjustment.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
v2 = zeros(14,3);
v2(1,1)=2;
v2(2,1)=2;
v2(3,1)=2;
v2(4,1)=2;
v2(5,1)=2;
v2(6,1)=2;
v2(7,1)=2;
v2(8,1)=5;
v2(9,1)=5;
v2(10,1)=5;
v2(11,1)=5;
v2(12,1)=6;
v2(13,1)=6;
v2(14,1)=6;
v2(1,2)=40;
v2(2,2)=46;
v2(3,2)=112;
v2(4,2)=47;
v2(5,2)=124;
v2(6,2)=119;
v2(7,2)=130;
v2(8,2)=27;
v2(9,2)=33;
v2(10,2)=99;
v2(11,2)=105;
v2(12,2)=1;
v2(13,2)=7;
v2(14,2)=73;
v2(1,3)=(-(y(11)*params(1)*y(10)*getPowerDeriv(y(4),params(1)-1,2)));
v2(2,3)=(-(y(11)*params(1)*T(5)));
v2(3,3)=v2(2,3);
v2(4,3)=(-(params(1)*y(10)*T(5)));
v2(5,3)=v2(4,3);
v2(6,3)=(-(params(1)*T(1)));
v2(7,3)=v2(6,3);
v2(8,3)=(-(T(2)*(T(7)*(-((-y(9))*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3))+T(6)*T(6)*T(8))));
v2(9,3)=(-(T(2)*(T(7)*(-1)/(y(3)*y(3))+T(6)*1/y(3)*T(8))));
v2(10,3)=v2(9,3);
v2(11,3)=(-(T(2)*1/y(3)*1/y(3)*T(8)));
v2(12,3)=(-(y(7)*getPowerDeriv(y(1),params(1),2)));
v2(13,3)=(-T(4));
v2(14,3)=v2(13,3);
g2 = sparse(v2(:,1),v2(:,2),v2(:,3),6,144);
end
