function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 13);

T = adjustment_nattinger_nogamma.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(8) = getPowerDeriv(y(1),params(1),1);
T(9) = getPowerDeriv(y(4),params(1)-1,1);
T(10) = params(1)*y(12)*T(9)+params(7)*((T(3)*(-params(3))*2*(y(11)-params(3)*y(4))-T(2)*2*2*y(4))/(T(3)*T(3))+(y(4)*params(3)*(-params(3))-params(3)*(y(11)-params(3)*y(4)))/(y(4)*y(4)));
T(11) = params(7)*(2*(y(11)-params(3)*y(4))/T(3)+params(3)/y(4));
T(12) = (-y(9))/(y(3)*y(3));
T(13) = getPowerDeriv(y(9)/y(3),(-params(6)),1);

end
