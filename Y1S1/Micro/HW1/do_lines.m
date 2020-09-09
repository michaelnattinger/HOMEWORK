function [lnp,itsct] = do_lines(pts1,pts2,prices)
% Calculates lines perp to prices going through the pts (takes pts1 as x
% and pts2 as y)
lnp = [0*pts1 0*pts2];
itsct = zeros(size(lnp,1)-1,size(lnp,2));
for i=1:length(pts1)
    lnp(i,1) = -(prices(i,1)/prices(i,2));
    lnp(i,2) = pts2(i) - lnp(i,1)*pts1(i);
end
for i=1:length(pts1)-1
    itsct(i,1) = (lnp(i+1,2) - lnp(i,2))/(lnp(i,1) - lnp(i+1,1));
    itsct(i,2) = lnp(i,1)*itsct(i,1) + lnp(i,2);
end