function [par, metaPar, txtPar] = pars_init_Oikopleura_dioica(metaData)

metaPar.model = 'std'; 

% reference parameter (not to be changed)
par.T_ref = C2K(20); free.T_ref = 0; units.T_ref = 'K';        label.T_ref = 'Reference temperature';

%% core primary parameters
par.z    = 0.0575*1.45;      free.z     = 1;   units.z     = '-';        label.z     = 'zoom factor';
par.F_m   = 4800;   free.F_m   = 0;   units.F_m   = 'l/d.cm^2'; label.F_m   = '{F_m}, max spec searching rate';
par.kap_X = 0.8;   free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve';
par.kap_P = 0.1;   free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces';
par.v     = 0.0223;  free.v     = 1;   units.v     = 'cm/d';     label.v     = 'energy conductance';
par.kap   = 0.06;   free.kap   = 1;   units.kap   = '-';        label.kap   = 'allocation fraction to soma';
par.kap_R = 0.95;  free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency';
par.p_M   = 1434;    free.p_M   = 1;   units.p_M   = 'J/d.cm^3'; label.p_M   = '[p_M], vol-spec somatic maint';
par.p_T   =  0;    free.p_T   = 0;   units.p_T   = 'J/d.cm^2'; label.p_T   = '{p_T}, surf-spec somatic maint';
par.k_J   = 0.001; free.k_J   = 0;   units.k_J   = '1/d';      label.k_J   = 'maturity maint rate coefficient';
par.E_G   = 2800*1.2;  free.E_G   = 1;   units.E_G   = 'J/cm^3';   label.E_G   = '[E_G], spec cost for structure';
par.E_Hb  = 2e-4*1.8;  free.E_Hb  = 1;   units.E_Hb  = 'J';        label.E_Hb  = 'maturity at birth';
par.E_Hp  = 1.3e-2*1.3;    free.E_Hp  = 1;   units.E_Hp  = 'J';        label.E_Hp  = 'maturity at puberty';
par.h_a   = 8e-7*80;  free.h_a   = 0;   units.h_a   = '1/d^2';    label.h_a   = 'Weibull aging acceleration';
par.s_G   = 15.53;  free.s_G   = 0;   units.s_G   = '-';        label.s_G   = 'Gompertz stress coefficient';

%% auxiliary parameters
par.T_A   = 11000;   free.T_A   = 0;    units.T_A = 'K';        label.T_A = 'Arrhenius temperature';
par.del_M = 0.266;   free.del_M = 1;    units.del_M = '-';      label.del_M = 'shape coefficient';

par.V_R = 0.0465;   free.V_R = 0;    units.V_R = '1/cm^3.d';      label.V_R = 'conversion factor for gonad length';

%% environmental parameters (temperatures are in auxData)
par.f = 1.0;        free.f     = 0;    units.f = '-';          label.f    = 'scaled functional response for 0-var data';
par.f1_tLT = 0.89;     free.f1_tLT  = 1;    units.f1_tLT = '-'; label.f_tLT = 'high food level Troes2002';
par.f2_tLT = 0.58;     free.f2_tLT  = 1;    units.f2_tLT = '-'; label.f_tLT = 'low food level Troes2002';

par.f1_Lomb2005 = 0.89;     free.f1_Lomb2005  = 1;    units.f1_Lomb2005 = '-'; label.f_Lomb2005 = 'high food level Lomb2005';
par.f2_Lomb2005 = 0.58;     free.f2_Lomb2005  = 1;    units.f2_Lomb2005 = '-'; label.f_Lomb2005 = 'low food level Lomb2005';


%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

par.d_V = 0.12;
par.d_E = 0.1679;
par.d_VC = 0.06;
par.d_EC = 0.084;

%% Pack output:
txtPar.units = units; txtPar.label = label; par.free = free; 

