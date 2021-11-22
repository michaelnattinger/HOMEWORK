function obj = smm_fm_3(x,MT,W,e,P)
P.prho = x(1);
P.psigma = x(2);
obj = smm_obj_3(MT,W,e,P);
end