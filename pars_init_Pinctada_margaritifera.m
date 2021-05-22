function [par, metaPar, txtPar] = pars_init_Pinctada_margaritifera(metaData)

metaPar.model = 'asj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 5523;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 0.30734;      free.z     = 0;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 57.73;    free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.78478;  free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.019731;     free.v     = 0;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.80981;    free.kap   = 0;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 8.3875;     free.p_M   = 0;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.0010789;  free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2370.151;   free.E_G   = 0;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 1.284e-04; free.E_Hb  = 0;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hs = 6.940e-04; free.E_Hs  = 0;   units.E_Hs = 'J';         label.E_Hs = 'maturity at settlement'; 
par.E_Hj = 2.615e+01; free.E_Hj  = 0;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 1.044e+03; free.E_Hp  = 0;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 5.784e-11;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.F_mb = 1.984;   free.F_mb  = 0;   units.F_mb = 'l/d.cm^2';  label.F_mb = '{F_mb}, max spec searching rate a larvae stage'; 
par.T_AH = 161000;    free.T_AH  = 0;   units.T_AH = 'K';         label.T_AH = 'Arrhenius temperature for upper boundary'; 
par.T_AL = 139.7;     free.T_AL  = 0;   units.T_AL = 'K';         label.T_AL = 'Arrhenius temperature for lower boundary'; 
par.T_H = 307.7;      free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'Upper temp boundary'; 
par.T_L = 286.2;      free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'Lower temp boundary'; 
par.del_M = 0.2781;   free.del_M = 0;   units.del_M = '-';        label.del_M = 'shape coefficient for shell heigth'; 
par.del_Mb = 1.0675;  free.del_Mb = 0;   units.del_Mb = '-';       label.del_Mb = 'shape coefficient shell diamaeter'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_TF = 0.2;   free.f_TF  = 0;   units.f_TF = '-';         label.f_TF = 'scaled functional response for TF data'; 
par.f_WdF = 0.2;  free.f_WdF = 0;   units.f_WdF = '-';     label.f_WdF = 'scaled functional response for WdF data';
par.f_WdJO = 0.1;     free.f_WdJO = 0;   units.f_WdJO = '-';       label.f_WdJO = 'scaled functional response for WdJO data'; 
par.f_WdJO1 = 0.9;    free.f_WdJO1 = 0;   units.f_WdJO1 = '-';      label.f_WdJO1 = 'scaled functional response for WdJO1 data'; 
par.f_lagoon = 0.4606;  free.f_lagoon = 0;   units.f_lagoon = '-';     label.f_lagoon = 'scaled functional response for S. Pouvreau in situ data'; 
par.f_tL = 0.31252;   free.f_tL  = 0;   units.f_tL = '-';         label.f_tL = 'scaled functional response for tL1 data'; 
par.f_tL2 = 1;   free.f_tL2 = 0;   units.f_tL2 = '-';        label.f_tL2 = 'scaled functional response for tL2 data, f = 1'; 
par.f_tL3 = 0.40308;  free.f_tL3 = 0;   units.f_tL3 = '-';        label.f_tL3 = 'scaled functional response for tL3 data'; 
par.f_tLC1 = 0.40405;  free.f_tLC1 = 0;   units.f_tLC1 = '-';       label.f_tLC1 = 'scaled functional response for tLC1 data'; 
par.f_tLC2 = 0.52221;  free.f_tLC2 = 0;   units.f_tLC2 = '-';       label.f_tLC2 = 'scaled functional response for tLC2 data'; 
par.f_tLC3 = 0.67177;  free.f_tLC3 = 0;   units.f_tLC3 = '-';       label.f_tLC3 = 'scaled functional response for tLC3 data'; 
par.f_tLC4 = 0.68758;  free.f_tLC4 = 0;   units.f_tLC4 = '-';       label.f_tLC4 = 'scaled functional response for tLC4 data'; 
par.p_Xm_amph = 1;    free.p_Xm_amph = 1;   units.p_Xm_amph = '';     label.p_Xm_amph = ''; 
par.p_Xm_isoc = 1;    free.p_Xm_isoc = 1;   units.p_Xm_isoc = '';     label.p_Xm_isoc = ''; 
par.p_Xm_nitz = 1;    free.p_Xm_nitz = 1;   units.p_Xm_nitz = '';     label.p_Xm_nitz = ''; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
