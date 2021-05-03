function [par, metaPar, txtPar] = pars_init_Pinctada_margaritifera(metaData)

metaPar.model = 'asj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 5523;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 0.13966;      free.z     = 0;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 55.41;   free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.041343;  free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.8782;   free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.0057272;    free.v     = 0;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9997;     free.kap   = 0;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 10.5033;    free.p_M   = 0;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.0014149;  free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2308.7102;  free.E_G   = 0;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 4.925e-09; free.E_Hb  = 0;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hs = 7.406e-08; free.E_Hs  = 0;   units.E_Hs = 'J';         label.E_Hs = 'maturity at settlement'; 
par.E_Hj = 1.268e-02; free.E_Hj  = 0;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 1.300e+00; free.E_Hp  = 0;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 5.784e-11;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.F_mb = 0.82677;   free.F_mb  = 1;   units.F_mb = 'l/d.cm^2';  label.F_mb = '{F_mb}, max spec searching rate a larvae stage'; 
par.F_m_amph = 0.3;   free.F_m_amph  = 1;   units.F_m_amph = 'l/d.cm^2';  label.F_m_amph = '{F_m.amph}, max spec searching rate with Amphidinium diet'; 
par.F_m_nitz = 0.3;   free.F_m_nitz  = 1;   units.F_m_nitz = 'l/d.cm^2';  label.F_m_nitz = '{F_m.nitz}, max spec searching rate with Nitzschia diet'; 
par.F_m_isoc = 0.3;   free.F_m_isoc  = 1;   units.F_m_isoc = 'l/d.cm^2';  label.F_m_isoc = '{F_m.isoc}, max spec searching rate with Isochrysis diet'; 
par.T_AH = 161000;  free.T_AH  = 0;   units.T_AH = 'K';         label.T_AH = 'Arrhenius temperature for upper boundary'; 
par.T_AL = 139.7;    free.T_AL  = 0;   units.T_AL = 'K';         label.T_AL = 'Arrhenius temperature for lower boundary'; 
par.T_H = 307.7;   free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'Upper temp boundary'; 
par.T_L = 286.2;     free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'Lower temp boundary'; 
par.del_M = 0.2781;  free.del_M = 0;   units.del_M = '-';        label.del_M = 'shape coefficient for shell heigth'; 
par.del_Mb = 1.678;   free.del_Mb = 0;   units.del_Mb = '-';       label.del_Mb = 'shape coefficient shell diamaeter'; 
par.f = 0.8;          free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
% par.f_WdJO = 0.3;     free.f_WdJO = 0;   units.f_WdJO = '-';       label.f_WdJO = 'scaled functional response for WdJO data'; 
% par.f_WdJO = 0.1;     free.f_WdJO = 0;   units.f_WdJO = '-';       label.f_WdJO = 'scaled functional response for WdJO data'; 
par.f_WdJO = 0.1;     free.f_WdJO = 0;   units.f_WdJO = '-';       label.f_WdJO = 'scaled functional response for WdJO data'; 
par.f_WdJO1 = 0.9;    free.f_WdJO1 = 0;   units.f_WdJO1 = '-';      label.f_WdJO1 = 'scaled functional response for WdJO1 data'; 
par.f_TF = 0.2;    free.f_TF = 0;   units.f_TF = '-';      label.f_TF = 'scaled functional response for TF data'; 
par.f_lagoon = 0.3;   free.f_lagoon = 1;   units.f_lagoon = '-';     label.f_lagoon = 'scaled functional response for S. Pouvreau in situ data'; 
par.f_tL = 0.3592;       free.f_tL  = 0;   units.f_tL = '-';         label.f_tL = 'scaled functional response for tL data'; 
par.f_tL2 = 0.8224;      free.f_tL2 = 0;   units.f_tL2 = '-';        label.f_tL2 = 'scaled functional response for tL2 data'; 
par.f_tLC1 = 0.21;    free.f_tLC1 = 0;   units.f_tLC1 = '-';       label.f_tLC1 = 'scaled functional response for tLC1 data'; 
par.f_tLC2 = 0.34;    free.f_tLC2 = 0;   units.f_tLC2 = '-';       label.f_tLC2 = 'scaled functional response for tLC2 data'; 
par.f_tLC3 = 0.76;    free.f_tLC3 = 0;   units.f_tLC3 = '-';       label.f_tLC3 = 'scaled functional response for tLC3 data'; 
par.f_tLC4 = 0.85;    free.f_tLC4 = 0;   units.f_tLC4 = '-';       label.f_tLC4 = 'scaled functional response for tLC4 data'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
