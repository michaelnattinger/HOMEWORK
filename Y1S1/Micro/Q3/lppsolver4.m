clear; close all; clc
% Solves for upper and lower bounds of quantities given a grid (4 by 4)
grid = [10  0  0  0; ... % enter grid here
        15 10  0  0; ...
        17 14 10  0; ...
        13 16 19 32];
    % solutions are in 'valid'
    % to check who matches with who,
    % ind`n' gives the column matched with row `n'
    
    
[~,n] = size(grid);
% solve for eqm - this is for 4d grid
tot = zeros(n,n,n,n);
for i = 1:n
    for  j= setdiff(1:n,i)
        for k=setdiff(1:n,[i j])
            for l = setdiff(1:n,[i j k])
                tot(i,j,k,l) = grid(1,i) + grid(2,j) + grid(3,k) + grid(4,l);
            end
        end
    end
end
eq = zeros(1,n);
[~,ind] = max(tot(:));
[ind1,ind2,ind3,ind4] = ind2sub(size(tot),ind); % optimal indices here
globmin = 0; %min(grid(:)) - 10; ZERO OUTSIDE OPTION
globmax = max(grid(:));
ops = globmin:globmax;
w = [0 0 0 0];
v = w; % V GOES TO LEFT, W GOES TO UP
valid = [];
for v1 = ops
    w(ind1) =  grid(1,ind1) - v1;
    for v2 = ops
        w(ind2) =  grid(2,ind2) - v2;
        for v3 = ops
            w(ind3) =  grid(3,ind3) - v3;
            for v4 = ops
                w(ind4) = grid(4,ind4) - v4;
            %w1 = w(1); w2 = w(2); w3 = w(3);
                v = [v1 v2 v3 v4];
                if min(w)>-1
                if check_lpp3(v,w,grid)
                    valid = [valid; v w];
                end
                end
            end
        end
    end
end
% we now have the full solution set.
