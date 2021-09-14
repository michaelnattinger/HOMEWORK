mkdir('pings')
clear; close all; clc
% PROGRAM NAME: 387vfigrowth.M
% This program generates the value function and decision rules for
% a nonstochastic growth model.
% Date: 2/17/03

% PARAMETERS
b=.99; %discount factor 
d=0.025; %depreciation rate
a=.36; %capital share
Z_g = 1.25; % technology
Z_b = 0.2;
P_gg = 0.977; % Transition probabilities
P_bb = 0.926;
P_gb = 0.023;
P_bg = 0.074;

% ASSET VECTOR % none of this changes
klb=0.01; %lower bound of grid points
%inc=0.025; %increments
kub=75;%upper bound of grid points
%k=[klb:inc:kub];% asset (row) vector
nkg = 3000; %<-- change this to 1000 for speed comparison w/ Julia/Fortran
k=linspace(klb,kub,nkg);
N=size(k);
N=N(1,2);
c1=ones(N,1); % column vector of ones
K=c1*k; % rows are k, cols are k'

% TABULATE CURRENT RETURN (UTILITY) FUNCTION
%% CHANGES HERE 9/9/21 Nattinger
% cs=zeros(N,N); %rows are k, cols are k'
% cs=K'.^a-(K-K'*(1-d));%cons
% is=find(cs<0);
% cs(is)=0;
% us=log(cs);
% t=isinf(us);
% j=find(t==1);
% us(j)=-realmax;
cs=zeros(N,N,2); %rows are k, cols are k'
cs(:,:,1)=Z_g*K'.^a-(K-K'*(1-d));%cons
cs(:,:,2)=Z_b*K'.^a-(K-K'*(1-d));

is=find(cs<0);
cs(is)=0;
us=log(cs);
t=isinf(us);
j=find(t==1);
us(j)=-realmax;

%%

% TABULATE INITIAL VALUE FUNCTION GUESS
visr=squeeze(us(:,1,:))'; %choose utility associated with k'=0 

pcntol=1; %tolerance for value function iteration
n=1; %if want to run vfi for a set number of iterations
tic
while pcntol >.0001
    %% MORE CHANGES IN HERE
   vis_g=c1*visr(1,:); %generates future value function matrix from above row vector
   vis_b = c1*visr(2,:);
   %CONSTRUCT TOTAL RETURN FUNCTION
%    wis=us+b*vis;
    wis_g=us(:,:,1)+b*(vis_g*P_gg + vis_b*P_gb);
    wis_b=us(:,:,2)+b*(vis_g*P_bg + vis_b*P_bb);
   %CHOOSE HIGHEST VALUE (ASSOCIATED WITH k' CHOICE)
   [vsr_g,I_g]=max(wis_g'); %since max gives highest element in each column of a matrix
   [vsr_b,I_b]=max(wis_b');
   
   n=n+1;
   vsr = [vsr_g; vsr_b];
   pcntol=max(max(abs(vsr-visr))); %use sup norm for tolerance
   %pcntol=tol/abs(vsr(1,N,1));
   visr=vsr; %update value functions
end
toc % display runtime in consol
% save 387vdr vsr I k;
% save 387parm b a d N inc klb kub;

%plot(k,k([I_g])-k)% plot change in the decision rule
%plot(k,vsr(1,:)) % plot value function

cd('pings')
% Make pretty graphs
figure
plot(k,vsr(1,:),'k')
hold on
plot(k,vsr(2,:),'r')
hold off
set(gcf,'Color',[1 1 1])
title('Value functions')
legend('Z^g','Z^b','Location','SouthEast')
saveas(gcf,'value.png')

figure
plot(k,k(I_g),'k')
hold on
plot(k,k(I_b),'r')
plot(k,k,'b')
hold off
set(gcf,'Color',[1 1 1])
title('Policy functions (K prime)')
legend('Z^g','Z^b','45 degree line','Location','SouthEast')
saveas(gcf,'policy.png')

figure
plot(k,k(I_g)-k,'k')
hold on
plot(k,k(I_b)-k,'r')
plot(k,0*k,'b')
hold off
set(gcf,'Color',[1 1 1])
title('Savings (K prime - K)')
legend('Z^g','Z^b','Location','SouthEast')
saveas(gcf,'savings.png')
cd('..')
