function [prdData, info] = predict_Oikopleura_dioica(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  % customized filters for allowable parameters of the standard DEB model (std)
  % for other models consult the appropriate filter function.
  filterChecks = k * v_Hp >= f1_tLT^3 || ...         % constraint required for reaching puberty with f1_tLT
                 ~reach_birth(g, k, v_Hb, f1_tLT) || ... % % constraint required for reaching birth with f1_tLT
             k * v_Hp >= f2_tLT^3 || ...        
                 ~reach_birth(g, k, v_Hb, f2_tLT) || ... 
             k * v_Hp >= f1_Lomb2005 || ...        
                 ~reach_birth(g, k, v_Hb, f1_Lomb2005) || ... 
             k * v_Hp >= f2_Lomb2005 || ...        
                 ~reach_birth(g, k, v_Hb, f2_Lomb2005);   
  
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_Ni = tempcorr(temp.Ni, T_ref, T_A);
  TC_tL_15 = tempcorr(temp.tL_S15, T_ref, T_A);
  TC_tL_20 = tempcorr(temp.tLtot_S20, T_ref, T_A);
  TC_XpA = tempcorr(temp.XpA, T_ref, T_A);
  TC_tJO = tempcorr(temp.tJO_f1, T_ref, T_A);
   
  % zero-variate data

  % life cycle
  pars_tp = [g k l_T v_Hb v_Hp];
  [t_p, t_b, l_p, l_b, info] = get_tp(pars_tp, f);
  
%   % initial
%   pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
%   U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M;                % cm, physical length at birth at f
  Wd_b = d_V * L_b^3 * (1 + f * w);       % g, dry weight at birth at f 
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth at f and T

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, physical length at puberty at f
  Wd_p = d_V * L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  aT_p = t_p/ k_M/ TC_ap;           % d, age at puberty at f and T

  % ultimate
  l_i = f - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate physical length at f
  Wd_i = d_V * L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % cum. reproduction at ultimate size at death
  pars_R = [kap; kap_R; g; k_J * TC_Ni; k_M  * TC_Ni; L_T; v  * TC_Ni; U_Hb/ TC_Ni; U_Hp/ TC_Ni]; % compose parameter vector at T
  time = tdeath.Ni;
  NT_i = cum_reprod(time, f, pars_R, L_b);
%   RT_i = TC_Ri * reprod_rate(L_i, f, pars_R);             % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
%% pack to output
  prdData.ab = aT_b;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wdb = Wd_b;
  prdData.Wdp = Wd_p;
  prdData.Wdi = Wd_i;
  prdData.Ni = NT_i;
  
%% uni-variate data
  
% time-length 15 degrees, troe2002
kT_J   =  k_J * TC_tL_15; kT_M = k_M  * TC_tL_15; vT = v * TC_tL_15;
UT_Hb  = U_Hb/ TC_tL_15; UT_Hp = U_Hp/ TC_tL_15; 
pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T
VT_R = V_R * TC_tL_15; % 1/cm^3.d, conversion factor for gonad length

% high food level
  [~, L] = cum_reprod(tL_S15(:,1), f1_tLT, pars_R);
  EL_S15 = L/del_M;
  [N, L, UT_E0] = cum_reprod(tLR_S15(:,1), f1_tLT, pars_R);
  ELR_S15 = VT_R * UT_E0 * N ./ (L .^ 2);

% low food level
  [~, L] = cum_reprod(tL_L15(:,1), f2_tLT, pars_R);
  EL_L15 = L/del_M;
  [N, L, UT_E0] = cum_reprod(tLR_L15(:,1), f2_tLT, pars_R);
  ELR_L15 = VT_R * UT_E0 * N ./ (L .^ 2);    
    
% 20 deg, total body length, two food levels
  % 15 degrees, troe2002
kT_J   =  k_J * TC_tL_20; kT_M = k_M  * TC_tL_20; vT = v * TC_tL_20;
UT_Hb  = U_Hb/ TC_tL_20; UT_Hp = U_Hp/ TC_tL_20; 
pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T
VT_R = V_R * TC_tL_20; % 1/cm^3.d, conversion factor for gonad length

[N, L, UT_E0] = cum_reprod(tLtot_S20(:,1), f1_tLT, pars_R);
LR = VT_R * UT_E0 * N ./ (L .^ 2);    
ELtot_S20 = LR + L/ del_M;  
  
[N, L, UT_E0] = cum_reprod(tLtot_L20(:,1), f2_tLT, pars_R);
LR = VT_R * UT_E0 * N ./ (L .^ 2);    
ELtot_L20 = LR + L/ del_M;  
  
% assimilation data
L = Lw.XpA * del_M; % cm, structural length
K = p_Am/ F_m; % J/d, half saturation coefficient 
F = XpA(:,1)./ (K + XpA(:,1)); % -, scaled functional response
EpA = F .* p_Am * L^2 * TC_XpA ; % J/d, assimilation rate

%% Lomb2005
% respiration and growth at two food levels
kT_J   =  k_J * TC_tJO ; kT_M = k_M  * TC_tJO ; vT = v * TC_tJO;
UT_Hb  = U_Hb/ TC_tJO; UT_Hp = U_Hp/ TC_tJO;
VT_R = V_R * TC_tJO; % 1/cm^3.d, conversion factor for gonad length
p_ref = p_Am * L_m^2; % max assimilation power at max size
X_gas = 1/ 24.4;     % M, mol of gas per litre at 20 C and 1 bar 
pars_R   = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T
pars_pow = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T_ref

% high food level
[N, L, UT_E0, L_b, L_p] = cum_reprod(tL_f1(:,1), f1_Lomb2005, pars_R);
l_b = L_b/ L_m; l_p = L_p/ L_m; % scale lengths at birth and puberty
  pACSJGRD = p_ref * scaled_power(L, f1_Lomb2005, pars_pow, l_b, l_p); % powers
  pADG = pACSJGRD(:,[1 7 5])';      % assimilation, dissipation, growth power
%   pADG(1,:) = 0; % exclude contributions from assimulation
  J_O = eta_O * pADG;               % J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M = - (n_M\n_O) * J_O;          % J_C, J_H, J_O, J_N in rows, A, D, G in cols
EJO_f1 = - TC_tJO * 1e6 * J_M(3,:)'/ 24/ X_gas;  % mul/h, dioxygen use
ELR_f1 = VT_R * UT_E0 * N ./ (L .^ 2) * 1e4;    % mum, gonad length
EL_f1 = L/ del_M * 1e4; % mum, trunk length

% low food level
[N, L, UT_E0, L_b, L_p] = cum_reprod(tL_f2(:,1), f2_Lomb2005, pars_R);
l_b = L_b/ L_m; l_p = L_p/ L_m; % scale lengths at birth and puberty
  pACSJGRD = p_ref * scaled_power(L, f2_Lomb2005, pars_pow, l_b, l_p); % powers
  pADG = pACSJGRD(:,[1 7 5])';      % assimilation, dissipation, growth power
%   pADG(1,:) = 0; % exclude contributions from assimulation
  J_O = eta_O * pADG;               % J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M = - (n_M\n_O) * J_O;          % J_C, J_H, J_O, J_N in rows, A, D, G in cols
EJO_f2 = - TC_tJO * 1e6 * J_M(3,:)'/ 24/ X_gas;  % mul/h, dioxygen use
ELR_f2 = VT_R * UT_E0 * N ./ (L .^ 2) * 1e4;    % mum, gonad length
EL_f2 = L/ del_M * 1e4; % mum, trunk length

%% pack to output
  prdData.tL_S15 = EL_S15;
  prdData.tLR_S15 = ELR_S15;
  prdData.tL_L15 = EL_L15;
  prdData.tLR_L15 = ELR_L15;
  prdData.tLtot_S20 = ELtot_S20;
  prdData.tLtot_L20 = ELtot_L20;
  prdData.XpA = EpA;
  prdData.tL_f1 = EL_f1;
  prdData.tLR_f1 = ELR_f1;
  prdData.tL_f2 = EL_f2;
  prdData.tLR_f2 = ELR_f2;
  prdData.tJO_f1 = EJO_f1;
  prdData.tJO_f2 = EJO_f2;
  
  
