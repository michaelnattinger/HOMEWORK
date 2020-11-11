function saddle_path = calc_saddle_p3(kgrid,kst,num,kss,Iss,del,R)

ng = length(kgrid);
grid = linspace(0,kst,num);
opt = 0*kgrid;
for i=1:ng
        k0 = grid*0 + kgrid(i);
        I0 = grid;
        for tt = 1:50
            k1 = I0 + (1-del)*k0;
            I1 = (R*I0 + k1 - kst)./(1-del);
            k0 = k1;
            I0 = I1;
        end
        [~,opt(i)] = min(abs(k1 - kss)+abs(I1 - Iss));
end
saddle_path = grid(opt);