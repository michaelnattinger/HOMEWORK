% PROGRAM NAME: ps1
%
% This program generates and plots the price dynamics given the first-order 
% difference equation discussed in the problem set given some initial price 
%
% Prepared by Fu Tan
% Last update: 09/08/15

clear;
clc;

%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%
r = 0.01; % interest rate
d = 1; % constant dividend
p0 = 100; % initial price
dim = 99; % terminal period t = 99

%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%
pvector = zeros(dim+1,1); % creating a vector of price from t=0 to t=99
tvector = linspace(0,dim,dim+1)'; % creating a vector for time from 0 to 99
pvector(1) = p0; % giving value to the first element of the price vector
                 % with the initial price

%%%%%%%%%%%
% DYNAMICS
%%%%%%%%%%%

for n = 2:dim+1 % starting from t = 1 to t = 100 
    pvector(n) = (1+r)*pvector(n-1)-d; % updating the price in the next period 
                                       % with the first-order difference
                                       % equation
end

%%%%%%%%
% PLOTS
%%%%%%%%
figure();
plot(tvector(:),pvector(:));
title('Price Dynamics');
xlabel('Time t'); ylabel('Price P_t');
legend('P_0 = 100','Location','Northeast');
axis([0 dim 0 150])
