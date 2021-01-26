clear
clc

% Solution code of Problem Set 2, ECON712
% Written by Anson Zhou, updated for 2017

%% Problem 1
z = 1;
alpha = 0.3;
delta = 0.1;
beta = 0.97;

% steady state
k_bar = ((1./(alpha*z))*(1./beta-(1-delta))).^(1./(alpha-1));
c_bar = z*k_bar.^alpha-delta*k_bar;

% compute the Jacobian matrix at steady state
J_11 = z*alpha*(k_bar).^(alpha-1)+(1-delta);
J_12 = -1;
J_21 = (c_bar*beta)*z*alpha*(alpha-1)*k_bar.^(alpha-2)...
    *(alpha*z*k_bar.^(alpha-1)+(1-delta));
J_22 = beta*(1-delta+z*alpha*k_bar^(alpha-1))...
    +(c_bar*beta)*z*alpha*(alpha-1)*k_bar.^(alpha-2)*(-1);
J = [J_11, J_12; J_21, J_22];

% eigenvalues and eigenvectors
[V,D]=eig(J);

% new steady state
z = z + 0.1;
k_bar2 = ((1./(alpha*z))*(1./beta-(1-delta))).^(1./(alpha-1));
c_bar2 = z*k_bar2.^alpha-delta*k_bar2;

% Jacobian matrix at new steady state
J2_11 = z*alpha*(k_bar2).^(alpha-1)+(1-delta);
J2_12 = -1;
J2_21 = (c_bar2*beta)*z*alpha*(alpha-1)*k_bar2.^(alpha-2)...
    *(alpha*z*k_bar2.^(alpha-1)+(1-delta));
J2_22 = beta*(1-delta+z*alpha*k_bar2^(alpha-1))...
    +(c_bar2*beta)*z*alpha*(alpha-1)*k_bar2.^(alpha-2)*(-1);
J2 = [J2_11, J2_12; J2_21, J2_22];

% eigenvalues and eigenvectors
[V2,D2]=eig(J2);

% determine coefficients for specific solution
t0 = 5;
m2 = (k_bar-k_bar2)/(V2(1,2).* D2(2,2).^t0);

% response at t0
c_t0 = c_bar2 + V2(2,2)*m2*(D2(2,2).^t0);

% trajectory using linearized saddle path
Prd = 20-t0+1;
Trj1 = zeros(2,Prd);
Trj1(1,1) = k_bar;
Trj1(2,1) = c_t0;
for i = 2:Prd
    Trj1(:,i)=[k_bar2; c_bar2]+J2*(Trj1(:,i-1)-[k_bar2; c_bar2]);
end

% trajectory using original equation and the "jump" of c_5 calculated under
% linearization
Trj2 = zeros(2,Prd);
Trj2(1,1)=k_bar;
Trj2(2,1)=c_t0;
for i = 2:Prd
    k = Trj2(1,i-1);
    c = Trj2(2,i-1);
    Trj2(1,i) = z*k.^alpha + (1-delta)*k - c;
    Trj2(2,i) = (1-delta+z*alpha*Trj2(1,i).^(alpha-1))*beta*c;
end

before = repmat([k_bar; c_bar],1, t0-1);
Trj1_whole = [before, Trj1];   % Transition under linear hypothesis
Trj2_whole = [before, Trj2];   % Actual transition using "jump" from linearization

time = 1:20;

% Linear plot
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj1_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj1_whole(2,:),'-*')
xlabel('Time')
title('c_t')


% Showing linear value under actual dynamics does not converge
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj2_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj2_whole(2,:),'-*')
xlabel('Time')
title('c_t')

% find actual c_t0 needed.
M = 10000;
C_range = linspace(c_bar, c_bar2, M);
Distance = zeros(1,M);
for i = 1:M
    Distance(1,i)=calib( C_range(i), alpha, beta, delta, z, k_bar, k_bar2, c_bar2);
end
[DD,I] = min(Distance);
c_star = C_range(I);

Trj3 = zeros(2,Prd);
Trj3(1,1)=k_bar;
Trj3(2,1)=c_star;
for i = 2:Prd
    k = Trj3(1,i-1);
    c = Trj3(2,i-1);
    Trj3(1,i) = z*k.^alpha + (1-delta)*k - c;
    Trj3(2,i) = (1-delta+z*alpha*Trj3(1,i).^(alpha-1))*beta*c;
end

Trj3_whole = [before, Trj3];  

% Actual transition using "shooting algorithm"
% Actual transition, using c_5 calculated from "shooting algorithm"
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(time,Trj3_whole(1,:),'-*')
xlabel('Time')
title('k_t')

subplot(2,1,2)       % add first plot in 2 x 1 grid
plot(time,Trj3_whole(2,:),'-*')
xlabel('Time')
title('c_t')

