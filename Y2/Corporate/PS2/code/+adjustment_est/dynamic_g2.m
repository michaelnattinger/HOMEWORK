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
    T = adjustment_est.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
v2 = zeros(27,3);
v2(1,1)=1;
v2(2,1)=1;
v2(3,1)=1;
v2(4,1)=2;
v2(5,1)=2;
v2(6,1)=2;
v2(7,1)=2;
v2(8,1)=2;
v2(9,1)=2;
v2(10,1)=2;
v2(11,1)=2;
v2(12,1)=2;
v2(13,1)=2;
v2(14,1)=2;
v2(15,1)=2;
v2(16,1)=2;
v2(17,1)=2;
v2(18,1)=5;
v2(19,1)=5;
v2(20,1)=5;
v2(21,1)=5;
v2(22,1)=6;
v2(23,1)=6;
v2(24,1)=6;
v2(25,1)=6;
v2(26,1)=6;
v2(27,1)=6;
v2(1,2)=1;
v2(2,2)=6;
v2(3,2)=81;
v2(4,2)=52;
v2(5,2)=60;
v2(6,2)=180;
v2(7,2)=61;
v2(8,2)=196;
v2(9,2)=62;
v2(10,2)=212;
v2(11,2)=174;
v2(12,2)=219;
v2(13,2)=188;
v2(14,2)=190;
v2(15,2)=220;
v2(16,2)=206;
v2(17,2)=221;
v2(18,2)=35;
v2(19,2)=41;
v2(20,2)=131;
v2(21,2)=137;
v2(22,2)=1;
v2(23,2)=6;
v2(24,2)=81;
v2(25,2)=7;
v2(26,2)=97;
v2(27,2)=86;
v2(1,3)=(-((-((y(1)*params(6)*(-params(3))-params(6)*(y(6)-params(3)*y(1)))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))));
v2(2,3)=(-((-params(6))/(y(1)*y(1))));
v2(3,3)=v2(2,3);
v2(4,3)=(-(y(14)*(params(1)*y(13)*getPowerDeriv(y(4),params(1)-1,2)+params(6)*((T(3)*T(3)*((-params(3))*2*(y(12)-params(3)*y(4))*2*2*y(4)+T(3)*(-params(3))*2*(-params(3))-((-params(3))*2*(y(12)-params(3)*y(4))*2*2*y(4)+T(2)*4))-(T(3)*(-params(3))*2*(y(12)-params(3)*y(4))-T(2)*2*2*y(4))*(T(3)*2*2*y(4)+T(3)*2*2*y(4)))/(T(3)*T(3)*T(3)*T(3))+(-((y(4)*params(3)*(-params(3))-params(3)*(y(12)-params(3)*y(4)))*(y(4)+y(4))))/(y(4)*y(4)*y(4)*y(4))))));
v2(5,3)=(-(y(14)*params(6)*((T(3)*2*(-params(3))-2*(y(12)-params(3)*y(4))*2*2*y(4))/(T(3)*T(3))+(-params(3))/(y(4)*y(4)))));
v2(6,3)=v2(5,3);
v2(7,3)=(-(y(14)*params(1)*T(9)));
v2(8,3)=v2(7,3);
v2(9,3)=(-T(10));
v2(10,3)=v2(9,3);
v2(11,3)=(-(1-params(3)));
v2(12,3)=v2(11,3);
v2(13,3)=(-(y(14)*params(6)*2/T(3)));
v2(14,3)=(-T(11));
v2(15,3)=v2(14,3);
v2(16,3)=(-(params(1)*T(1)));
v2(17,3)=v2(16,3);
v2(18,3)=(-(T(5)*(T(13)*(-((-y(9))*(y(3)+y(3))))/(y(3)*y(3)*y(3)*y(3))+T(12)*T(12)*T(14))));
v2(19,3)=(-(T(5)*(T(13)*(-1)/(y(3)*y(3))+T(12)*1/y(3)*T(14))));
v2(20,3)=v2(19,3);
v2(21,3)=(-(T(5)*1/y(3)*1/y(3)*T(14)));
v2(22,3)=(-(y(7)*getPowerDeriv(y(1),params(1),2)-(2*y(1)*2*y(1)*2*y(1)*params(6)*(-params(3))*2*(-params(3))-(2*y(1)*params(6)*(-params(3))*2*(y(6)-params(3)*y(1))-2*T(7))*(2*2*y(1)+2*2*y(1)))/(2*y(1)*2*y(1)*2*y(1)*2*y(1))));
v2(23,3)=(2*y(1)*params(6)*2*(-params(3))-2*params(6)*2*(y(6)-params(3)*y(1)))/(2*y(1)*2*y(1));
v2(24,3)=v2(23,3);
v2(25,3)=(-T(8));
v2(26,3)=v2(25,3);
v2(27,3)=2*params(6)/(2*y(1));
g2 = sparse(v2(:,1),v2(:,2),v2(:,3),7,256);
end
