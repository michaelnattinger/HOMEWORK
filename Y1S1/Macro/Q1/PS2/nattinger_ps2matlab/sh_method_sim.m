ng = 50000;
grid = linspace(1.099772720277396,ss.y.c,ng);
y0 = [y0(1)+0*grid;grid];
T = 50;
exact_sim
[~,i_opt] = min(abs(ss.y.c - simsim(2,end,:)));
sh_sim = squeeze(simsim(:,:,i_opt));
sh_c = grid(i_opt);