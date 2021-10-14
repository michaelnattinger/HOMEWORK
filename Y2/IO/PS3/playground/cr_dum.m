function f=cr_dum(long_id)

% cr_dum	This function creates a set of dummies for each of the values defined by long_id

b = sort(long_id);
b1 = [1;diff(b)];
b2 = b(b1>0);
clear b1 b
f = sparse(zeros(size(long_id,1), size(b2,1)));
for i = 1:size(b2,1)
	f(:,i) = sparse((long_id==b2(i)));
end
