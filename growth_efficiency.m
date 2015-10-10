% tu peux utiliser ce code pour voir plus facilement quel est to kap_G
% selon les valeur que tu mets dans pars_init:


[data, auxData, metaData, txtData, weights] = mydata_Oikopleura_dioica;
[par, metaPar, txtPar] = pars_init_Oikopleura_dioica(metaData);
 
 
 % kap_G = mu_V * d_V/ w_V/ E_G; -, growth conversion efficiency
 
 
 cp = parscomp_st(par); % computes compound parameters;
 
 fprintf('growth conversion efficiency kap_G is %1.2f \n',cp.kap_G) % print value of growth conversion efficiency to screen