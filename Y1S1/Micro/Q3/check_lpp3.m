function res = check_lpp3(v,w,grid)
res = 1;
n = length(v);
for i=1:n
    for j=1:n
        if grid(i,j)>v(i)+w(j); res = 0; return; end 
    end
end
end