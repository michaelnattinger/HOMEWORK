function obj = smm_fm_2(x,MT,W,e,P)
P.prho = x(1);
P.psigma = x(2);
obj = smm_obj_2(MT,W,e,P);
end