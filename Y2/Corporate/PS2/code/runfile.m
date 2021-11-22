clear; close all; clc
%% First part: run (and plot) impulse responses
% run dynare
dynare adjustment_nattinger.mod

% collect results into matrix irfs, normalize to percent deviation, and plot
irfs = [k_eps_z q_eps_z i_eps_z z_eps_z sdf_eps_z c_eps_z];
vars = oo_.var_list;
ssvals = oo_.steady_state';
irfs = irfs./ssvals;
sims = oo_.endo_simul;
xxx = 1:40;

figure
for ii=1:6
    subplot(3,2,ii)
    plot(xxx,100*irfs(:,ii),'b--')
    hold on
    plot(xxx,0*xxx,'k')
    hold off
    title(['response of ' vars{ii} ' to technology shock'])
    ylabel('percent')
end
set(gcf,'Color',[1 1 1])
suptitle('Impulse responses, \gamma = 2')
cd('..\pings')
saveas(gcf,'fig1.png')
cd('..\code')

% Now we do the same thing with a different gamma value
dynare adjustment_nattinger_nogamma.mod
% collect results into matrix irfs, normalize to percent deviation, and plot
irfs = [k_eps_z q_eps_z i_eps_z z_eps_z sdf_eps_z c_eps_z];
vars = oo_.var_list;
ssvals = oo_.steady_state';
irfs = irfs./ssvals;

figure % wow do these look different. 10% drop in consumption on impact to
       % invest in capital and take advantage of high z values.
for ii=1:6
    subplot(3,2,ii)
    plot(xxx,100*irfs(:,ii),'b--')
    hold on
    plot(xxx,0*xxx,'k')
    hold off
    title(['response of ' vars{ii} ' to technology shock'])
    ylabel('percent')
end
set(gcf,'Color',[1 1 1])
suptitle('Impulse responses, \gamma = 0')
cd('..\pings')
saveas(gcf,'fig2.png')
cd('..\code')

T = size(sims,1);
YY = sims(3:end,strcmp(vars,'i'))./sims(2:end-1,strcmp(vars,'k'));
XX = [ones(T-2,1) sims(2:end-1,strcmp(vars,'q'))  (sims(2:end-1,strcmp(vars,'z')).* sims(1:end-2,strcmp(vars,'k')).^(0.7))./sims(1:end-2,strcmp(vars,'k'))]; 
BB = XX\YY; % ols coefficients

vars = {'alpha' 'beta 1' 'beta 2'}';
tab = table(vars,BB);
cd('..')
table2latex_mid(tab,'ols_table')
cd('code')


%% Part 2: estimation

KK = sims(:,strcmp(vars,'k'));
QQ = sims(:,strcmp(vars,'q'));
%KK_obs = KK; % no measurement error here
QQ_obs = QQ + 0.01*normrnd(0,1,T,1); % measurement error here
k = KK;
qob = QQ_obs;
save('data_mat.mat','k','qob')

% use dynare to estimate parameters
dynare adjustment_est.mod
dynare adjustment_est_D.mod








