clear; close all; clc
n = 30;                 % # individuals
buyers = 2:2:n;         % which are buyers?
sellers = 1:2:n;        % which are sellers?
v_buyers = 2*buyers;    % buyer valuations
v_sellers = sellers.^2;  % seller valuations
q_buyers = 2+0*buyers;  % quantity each buyer wants to buy
q_sellers = 3+0*buyers; % quantity each supplier wants to sell
% returns market clearing prices in mktclr % note: does not return edge
% cases, so most useful to get quantities

pmax = max([v_buyers v_sellers]);
pmin = min([v_buyers v_sellers]);
p = 0:0.5:pmax;
q = [];
mktclr = [];
for pn = p
    i_d = find(v_buyers>pn);
    i_s = find(v_sellers<pn);
%     ind_d = length(find(v_buyers==pn));
%     ind_s =length(find( v_sellers == pn));
%     if i_d <= 1; ind_d = 0; end
%     if i_s <= 1; ind_s = 0; end
%     if (length(i_d) == length(i_s))|| (length(i_d)+ind_d == length(i_s)) ||(length(i_d) == length(i_s))+ind_s
    if (length(i_d) == length(i_s))
        mktclr = [mktclr pn];
        q = [q length(i_d)];
    end
end