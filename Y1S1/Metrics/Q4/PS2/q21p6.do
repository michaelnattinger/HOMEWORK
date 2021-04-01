pause on
use "C:\Users\micha\OneDrive\Documents\HOMEWORK\Y1S1\Metrics\Q4\PS3\LM2007.dta", clear
gen T = (povrate60>=59.2)
gen Tnx = (povrate60 - 59.2)*T
* h = 8 is baseline
reg mort_age59_related_postHS povrate60 Tnx T if povrate60>=45.4&povrate60<=72,r 
reg mort_age59_related_postHS povrate60 Tnx T if povrate60>=52.3&povrate60<=66.1,r 
reg mort_age59_related_postHS povrate60 Tnx T if povrate60>=38.4&povrate60<=80,r

reg mort_age25plus_related_postHS povrate60 Tnx T if povrate60>=45.4&povrate60<=72,r 
reg mort_age25plus_related_postHS povrate60 Tnx T if povrate60>=52.3&povrate60<=66.1,r 
reg mort_age25plus_related_postHS povrate60 Tnx T if povrate60>=38.4&povrate60<=80,r