function obj = smm_fm(x,MT,W,e,P)
P.prho = x(1);
P.psigma = x(2);
obj = smm_obj(MT,W,e,P);
end