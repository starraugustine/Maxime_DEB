function [Edata, EL_S15, GL_S15,...
EL_L15, GL_L15,...
BL_S20,...
BL_L20,...
ASS,eJO,eJO_L] = predict_Oikopleura_dioica(par, Data, tL_S15, gL_S15,...
tL_L15, gL_L15,...
bL_S20,...
bL_L20,...
Ass,O2_H,O2_L)

% [Edata, EL] = predict_Oikopleura_dioica(par, data, tL)
% par: (19)-vector with parameters (see below)
% data:(r_d,1)-matrix with zero-variate data (not some elements are used)
% tL: (r_L,1)-matrix with uni-variate data with time since birth, length
% Edata:(r_d,1)-matrix with expected values for data(:,2)
% EL:(r_L,1)-matrix with expected values for tL(:,2)
  
  global n_O n_M dwm ass_par CONC X_gas L_O2 L_O2_L% import values set in mydata_Oikopleura_dioica

%% unpack par
  T_ref  = par(1); % K, temp for which rate pars are given 
  T_A    = par(2); % K, Arrhenius temp
  f      = par(3); % scaled functional response
  z      = par(4); % zoom factor
  del_M  = par(5); % -, shape coefficient to convert vol-length to physical length
  F_m    = par(6); % l/d.cm^2, {F_m} max spec searching rate
  kap_X  = par(7); % -, digestion efficiency of food to reserve
  v      = par(8); % cm/d, energy conductance
  kap    = par(9); % -, alloaction fraction to soma = growth + somatic maintenance
  kap_R  = par(10);% -, reproduction efficiency
  p_M    = par(11);% J/d.cm^3, [p_M] vol-specific somatic maintenance
  p_T    = par(12);% J/d.cm^2, {p_T} surface-specific som maintenance
  k_J    = par(13);% 1/d, < k_M = p_M/E_G, maturity maint rate coefficient
  E_G    = par(14);% J/cm^3, [E_G], spec cost for structure
  E_Hb   = par(15);% J, E_H^b maturity level at birth
  E_Hp   = par(16);% J, E_H^b maturity level at birth
  h_a    = par(17);% 1/d^2, Weibull aging acceleration
  s_G    = par(18);% -, Gompertz stress coefficient
  v_R    = par(19);% cm/d, contribution of gonads to trunc length
  kap_X_P= par(20);% Sdefecation efficiency

  d_V = dwm(2,1);  d_E = dwm(3,1);  d_X = dwm(1,1);  % g/cm^3, specific densities for structure and reserve
  w_V = dwm(2,2);  w_E = dwm(3,2);  w_X = dwm(1,2);  % g/mol,  molecular weights for structure and reserve
  mu_V = dwm(2,3); mu_E = dwm(3,3); mu_X = dwm(1,3); mu_P = dwm(4,3); % J/mol,  chemical potentials for structure and reserve
  d_VC = dwm(2,4); d_EC = dwm(3,4); d_XC = dwm(1,4); % g C /cm^3, specific densities for structure and reserve

  % Selected copy-paste from parscomp & statistics
  p_Am = z * p_M/ kap;             % J/d.cm^2, {p_Am} max spec assimilation flux
  J_E_Am = p_Am/ mu_E;             % mol/d.cm^2, {J_EAm}, max surface-spec assimilation flux
  k_M = p_M/ E_G;                  % 1/d, somatic maintenance rate coefficient
  k = k_J/ k_M;                    % -, maintenance ratio
  % p_Xm = p_Am/ kap_X;            % J/d.cm^2, max spec feeding power

  M_V = d_V/ w_V;                  % mol/cm^3, [M_V] volume-specific mass of structure
  kap_G = M_V * mu_V/ E_G;         % -, growth efficiency
  y_V_E = mu_E * M_V/ E_G;         % mol/mol, yield of structure on reserve
  y_E_V = 1/ y_V_E;                % mol/mol, yield of reserve on structure
 
 %  Respiration & Feeding processes
  p_Xm     = p_Am/ kap_X;          % J/d.cm^2, max spec feeding power
  J_X_Am   = p_Xm/ mu_X;           % mol/d.cm^2, {J_XAm}, max surface-spec ingestion rate
  y_E_X    = kap_X * mu_X/ mu_E;   % mol/mol, yield of reserve on food
  y_X_E    = 1/ y_E_X;             % mol/mol, yield of food on reserve
  K        = J_X_Am/ F_m;          % c-mol X/l, half-saturation coefficient
  y_P_X    = kap_X_P * mu_X/ mu_P;    % mol/mol, yield of faeces on food 
  y_P_E    = y_P_X/ y_E_X;            % mol/mol, yield of faeces on reserve
  y_V_E    = mu_E * M_V/ E_G;         % mol/mol, yield of structure on reserve
  y_E_V    = 1/ y_V_E;                % mol/mol, yield of reserve on structure
  
  E_m = p_Am/ v;                   % J/cm^3, [E_m] reserve capacity 
  m_Em = y_E_V * E_m/ E_G;         % mol/mol, reserve capacity 
  g = E_G/ kap/ E_m;               % -, energy investment ratio
  w = m_Em * w_E/ w_V;             % -, contribution of reserve to weight

  L_m = v/ k_M/ g;                 % cm, maximum structural length
  L_T = p_T/ p_M;                  % cm, heating length (also applies to osmotic work)
  l_T = L_T/ L_m;                  % -, scaled heating length
  L_i = (f - l_T) * L_m;           % cm, ultimate structural length

  % maturity at birth = puberty
  M_Hb = E_Hb/ mu_E;               % mmol, maturity at birth  
  U_Hb = M_Hb/ J_E_Am;             % cm^2 d, scaled maturity at birth
  u_Hb = U_Hb * g^2 * k_M^3/ v^2;  % -, scaled maturity at birth
  V_Hb = U_Hb/ (1 - kap);          % cm^2 d, scaled maturity at birth
  v_Hb = V_Hb * g^2 * k_M^3/ v^2;  % -, scaled maturity at birth

    % maturity at puberty
  M_Hp = E_Hp/ mu_E;               % mol, maturity at puberty
  U_Hp = M_Hp/ J_E_Am;             % cm^2 d, scaled maturity at puberty 
  u_Hp = U_Hp * g^2 * k_M^3/ v^2;  % -, scaled maturity at puberty  
  V_Hp = U_Hp/ (1 - kap);          % cm^2 d, scaled maturity at puberty
  v_Hp = V_Hp * g^2 * k_M^3/ v^2;  % -, scaled maturity at puberty
  
   %  Mass-power couplers
  eta_XA = y_X_E/mu_E;             % mol/J, food-assim energy coupler
  eta_PA = y_P_E/mu_E;             % mol/J, faeces-assim energy coupler
  eta_VG = y_V_E/mu_E;             % mol/J, struct-growth energy coupler
  eta_O = [-eta_XA      0             0;              % mol/J, mass-energy coupler
	          0           0        eta_VG;              % used in: J_O = eta_O * p
	          1/mu_E    -1/mu_E   -1/mu_E;
	          eta_PA     0              0];
  
  %% zero-variate data

  % birth
  pars_lb = [g; k; v_Hb];                 % compose parameter vector
  [t_b, l_b, info] = get_tb(pars_lb, f);  % -, scaled and length at birth at f
  if info ~= 1                            % numerical procedure failed
    fprintf('warning: invalid parameter value combination for get_tb \n')
  end
  
  L_b = L_m * l_b;                        % cm, structural length at birth at f
  Lw_b = L_b/ del_M; 
  W_b = 1e6 * L_b^3 * d_VC * (1 + f * w); % mug, dry weight at birth at f
  TC = tempcorr(Data(1,1), T_ref, T_A);   % -, Temperature Correction factor for a_b; ONLY WORK WITH T_A
  kT_M = k_M * TC;                        % 1/d, correct k_M for temperature
  aT_b = t_b/ kT_M;  

  % puberty and ultimate size
  TC = tempcorr(Data(2,1), T_ref, T_A);
  pars_tp = [g; k; l_T; v_Hb; v_Hp];                % compose parameter vector
  [t_p t_b l_p l_b info] = get_tp(pars_tp, f, l_b); % -, scaled length at birth at f
  if info ~= 1                                      % numerical procedure failed
      fprintf('warning: invalid parameter value combination for get_tp \n')
  end
  L_p = L_m * l_p;                                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;
  W_p = 1e6 * L_p^3 * d_VC * (1 + f * w);
  kT_M = k_M * TC;
  aT_p = t_p/ kT_M;                                 % d, age at puberty at f and T


  % length at death at t = 6.625 d -->Data(10,2)
  TC = tempcorr(Data(10,1), T_ref, T_A);   
  kT_M = k_M * TC; kT_J = k_J * TC;                   % 1/d, correct k_M and k_J for temperature
  vT = v * TC;                                        % cm/d, energy conductance at T
  UT_Hb = U_Hb/ TC;
  UT_Hp = U_Hp/ TC;
  ir_B = 3/ kT_M + 3 * f * L_m/ vT; r_B = 1/ ir_B;    % d, 1/von Bert growth rate
  L_8 = L_i - (L_i - L_b) * exp( - r_B * Data(10,2)); % cm, expected length at time
  Lw_8 = L_8/ del_M ;                                 % cm, total trunc length (see uni-var data)
 
 % weight and reproduction at death at t = 6.625 d -->Data(10,2)
  pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T
%  [N_8 UE0] = cum_reprod(Data(10,2), f, pars_R, L_b);               % #, cumulative number of eggs at 8 d, UE0 
[N_8 L_8 UE0] = cum_reprod(Data(10,2), f, pars_R, L_b);
% This assumes that reproduction overheads are paid at the conversion of the buffer to eggs
  W_8 = 1e6 * L_8^3 * d_VC * (1 + f * w);                      % mug, dry weight at f for t = data(10,1) 
  Lw_8 = L_8/ del_M ;                                          % cm, total trunc length (see uni-var data)


  GL_8 = v_R * UE0 * N_8 ./ L_8 .^ 2;
  
  % UEO-->egg scaled reserve, N_8-->egg produced at death, product= scaled total energy invested in reprod at death 
  %dim(UE0) = time * length^2; dim(v_R) = length/ time
  % life span; death is not by ageing, but non-the-less described by it
  pars_tm = [g; k; l_T; v_Hb; v_Hp; h_a/ k_M^2; s_G];      % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b, l_p);                    % -, scaled mean life span at T_ref
  TC = tempcorr(Data(10,1), T_ref, T_A);                   % -, Temperature Correction factor fpr a_m; ONLY WORK WITH T_A
  kT_M = k_M * TC;                                         % 1/d, correct k_M for temperature
  aT_m = t_m/ kT_M;                                        % d, mean life span at 
  % pack output for zero-variate data
  Edata = [aT_b; aT_p; Lw_b; Lw_p; Lw_8; W_b; W_p; W_8; N_8; aT_m; % real data, T-corrected
           v; kap; kap_R; p_M; p_T; k_J; kap_G]; % pseudo data at T_ref

  %%%%%%%%%%%%% uni-variate data: total trunc length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = 273 + 15; f = 0.89;%1;                                            % this information was given in mydata_my_pet
  TC = tempcorr(T, T_ref, T_A);                                         % -, Temperature Correction factor; ONLY WORK WITH T_A
  kT_M = k_M * TC; kT_J = k_J * TC;                                     % 1/d, correct k_M for temperature  
  vT = v * TC;                                                          % cm/d, correct v for temperature  
  UT_Hb = U_Hb/ TC;  d_V = dwm(2,1);  d_E = dwm(3,1);  d_X = dwm(1,1);  % g/cm^3, specific densities for structure and reserve
  w_V = dwm(2,2);  w_E = dwm(3,2);  w_X = dwm(1,2);                     % g/mol,  molecular weights for structure and reserve
  mu_V = dwm(2,3); mu_E = dwm(3,3); mu_X = dwm(1,3);                    % J/mol,  chemical potentials for structure and reserve
  d_VC = dwm(2,4); d_EC = dwm(3,4); d_XC = dwm(1,4);                    % g C /cm^3, specific densities for structure and reserve
  % Selected copy-paste from parscomp & statistics
  p_Am = z * p_M/ kap;             % J/d.cm^2, {p_Am} max spec assimilation flux
  J_E_Am = p_Am/ mu_E;             % mol/d.cm^2, {J_EAm}, max surface-spec assimilation flux
  k_M = p_M/ E_G;                  % 1/d, somatic maintenance rate coefficient
  k = k_J/ k_M;                    % -, maintenance ratio

  M_V = d_V/ w_V;                  % mol/cm^3, [M_V] volume-specific mass of structure
  kap_G = M_V * mu_V/ E_G;         % -, growth efficiency
  y_V_E = mu_E * M_V/ E_G;         % mol/mol, yield of structure on reserve
  y_E_V = 1/ y_V_E;                % mol/mol, yield of reserve on structure
 
 %  Respiration & Feeding processes
  p_Xm     = p_Am/ kap_X;                          % J/d.cm^2, max spec feeding power
  J_X_Am   = p_Xm/ mu_X;                           % mol/d.cm^2, {J_XAm}, max surface-spec ingestion rate
  y_E_X    = kap_X * mu_X/ mu_E;                   % mol/mol, yield of reserve on food
  y_X_E    = 1/ y_E_X;                             % mol/mol, yield of food on reserve
  K        = J_X_Am/ F_m;                          % c-mol X/l, half-saturation coefficient
  UT_Hp = U_Hp/ TC;
  ir_B = 3/ kT_M + 3 * f * L_m/ vT; r_B = 1/ ir_B; % d, 1/von Bert growth rate

  VT_Hb = UT_Hb/ (1 - kap);
  v_Hb = VT_Hb * g^2 * kT_M^3/ vT^2;
  pars_lb = [g; k; v_Hb];
  l_b = get_lb(pars_lb, f);                        % -, scaled length at birth at f
  L_b = l_b * L_m;  L_i = f * L_m - L_T;

  L = (L_i - (L_i - L_b) * exp( - r_B * tL_S15(:,1)));         % cm, structural length
  pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T_ref
%  [N UE0] = cum_reprod(tL_S15(:,1), f, pars_R, L_b);           % #, cumulative number of eggs, UE0 
[N L UE0] = cum_reprod(tL_S15(:,1), f, pars_R, L_b);

  EL_S15 = L/ del_M ;                                          % cm, total trunc length
  vT_R = v_R*TC;
  GL_S15 = vT_R * UE0 * N ./ (L .^ 2);
%  [N UE0] = cum_reprod(tL_S15(:,1), 1, pars_R, L_b);           % #, cumulative number of eggs, UE0 
[N L UE0] = cum_reprod(tL_S15(:,1), f, pars_R, L_b);

  E0  = UE0 * p_Am ;                                           % J, initial energy in egg
  L_0001 = L_i - (L_i - L_b) * exp( - r_B * 0.001);
  % the cumulative reproductive material is deposited on (a fixed fraction of) body surface
  % so it contributes to length proportional to N/ L^2
  % dim(UE0) = time * length^2; dim(v_R) = length/ time
  
  %%%%%%%% uni-variate data: total trunc length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = 273 + 15; f = 0.58;                          % this information was given in mydata_my_pet
  TC = tempcorr(T, T_ref, T_A);                    % -, Temperature Correction factor; ONLY WORK WITH T_A
  kT_M = k_M * TC; kT_J = k_J * TC;                % 1/d, correct k_M for temperature  
  vT = v * TC;                                     % cm/d, correct v for temperature  
  UT_Hb = U_Hb/ TC;
  UT_Hp = U_Hp/ TC;
  ir_B = 3/ kT_M + 3 * f * L_m/ vT; r_B = 1/ ir_B; % d, 1/von Bert growth rate

  VT_Hb = UT_Hb/ (1 - kap);
  v_Hb = VT_Hb * g^2 * kT_M^3/ vT^2;
  pars_lb = [g; k; v_Hb];
  l_b = get_lb(pars_lb, 1);                        % -, scaled length at birth at f
  L_b = l_b * L_m;  L_i = f * L_m - L_T;
  
  L = (L_i - (L_i - L_b) * exp( - r_B * tL_L15(:,1)));         % cm, structural length
  pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T_ref
%  [N UE0] = cum_reprod(tL_L15(:,1), f, pars_R, L_b);           % #, cumulative number of eggs, UE0 
  [N L UE0] = cum_reprod(tL_L15(:,1), f, pars_R, L_b);
  
  EL_L15 = L/ del_M ;                                          % cm, total trunc length
  vT_R = v_R*TC;
  GL_L15 = vT_R * UE0 * N ./ L .^ 2;
  
  %%%%%%%%%%% uni-variate data: total trunc length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = 273 + 20; f = 0.89;                             % this information was given in mydata_my_pet
  TC = tempcorr(T, T_ref, T_A);                       % -, Temperature Correction factor; ONLY WORK WITH T_A
  kT_M = k_M * TC; kT_J = k_J * TC;                   % 1/d, correct k_M for temperature  
  vT = v * TC;                                        % cm/d, correct v for temperature  
  UT_Hb = U_Hb/ TC;
  UT_Hp = U_Hp/ TC;
  ir_B = 3/ kT_M + 3 * f * L_m/ vT; r_B = 1/ ir_B;    % d, 1/von Bert growth rate

  VT_Hb = UT_Hb/ (1 - kap);
  v_Hb = VT_Hb * g^2 * kT_M^3/ vT^2;
  pars_lb = [g; k; v_Hb];
  l_b = get_lb(pars_lb, f);                           % -, scaled length at birth at f
  L_b = l_b * L_m;  L_i = f * L_m - L_T;
    
 L = (L_i - (L_i - L_b) * exp( - r_B * bL_S20(:,1)));          % cm, structural length
  pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T_ref
%  [N UE0] = cum_reprod(bL_S20(:,1), f, pars_R, L_b);
  [N L UE0] = cum_reprod(bL_S20(:,1), f, pars_R, L_b);

  EL = L/ del_M ;                                              % cm, total trunc length
  vT_R = v_R*TC;
  GL = vT_R *UE0 * N ./ L .^ 2;
  BL_S20 = EL+GL;
  
  %%%%%%%%%%% uni-variate data: total trunc length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = 273 + 20; f = 0.58;                          % this information was given in mydata_my_pet
  TC = tempcorr(T, T_ref, T_A);                    % -, Temperature Correction factor; ONLY WORK WITH T_A
  kT_M = k_M * TC; kT_J = k_J * TC;                % 1/d, correct k_M for temperature  
  vT = v * TC;                                     % cm/d, correct v for temperature  
  UT_Hb = U_Hb/ TC;
  UT_Hp = U_Hp/ TC;
  ir_B = 3/ kT_M + 3 * f * L_m/ vT; r_B = 1/ ir_B; % d, 1/von Bert growth rate

  VT_Hb = UT_Hb/ (1 - kap);
  v_Hb = VT_Hb * g^2 * kT_M^3/ vT^2;
  pars_lb = [g; k; v_Hb];
  l_b = get_lb(pars_lb, 1);                        % -, scaled length at birth at f
  L_b = l_b * L_m;  L_i = f * L_m - L_T;
    
  L = (L_i - (L_i - L_b) * exp( - r_B * bL_L20(:,1)));         % cm, structural length
  pars_R = [kap; kap_R; g; kT_J; kT_M; L_T; vT; UT_Hb; UT_Hp]; % compose parameter vector at T_ref
% [N UE0] = cum_reprod(bL_L20(:,1), f, pars_R, L_b);            % #, cumulative number of eggs, UE0 
 [N L UE0] = cum_reprod(bL_L20(:,1), f, pars_R, L_b);
 
  EL = L/ del_M ;                                              % cm, total trunc length
  vT_R = v_R*TC;
  GL = vT_R * UE0 * N ./ L .^ 2;
  
  BL_L20 = EL+GL;
    
 %%%%%%%%%%%%%%%%%%%%%%%% Assimilation data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 T_f  = exp(((T_A)/(T_ref))-((T_A)/(ass_par(2)+273)));
 k_a  = K*mu_X ;                                               %Cmol/l * J/mol
 food =CONC*ass_par(1)./(CONC*ass_par(1)+k_a);   
 ASS = p_Am * T_f * (((ass_par(3))*del_M)^2).*food;
 
 %%%%%%%%%%%%%%%%%%%%%%%% Respiration rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calls subfunction *mineral_fluxes* f=1
f=1;
%L0       = (Wd0/conversion/d_V/ (1 + f * w)).^(1/ 3);  % cm, structural length
L0         = L_O2*del_M;
pars_tj    = [g; k; l_T; v_Hb; v_Hp]; 
temp_pars  = T_A;
% calculate mineral fluxes for an individual of structural length L0
 J_M  = mineral_fluxes(L0, f, L_m, E_m, p_M, E_G, kap, k_J, v, p_Am, p_T, E_Hp, kap_R, eta_O, n_M, n_O, pars_tj);
         
% dioxygen consumption
   TC    = tempcorr(15+273, T_ref, temp_pars);  % -, temperature correction factor
   eJO   = -1e6 * J_M(3,:)' .* TC* X_gas/ 24;           % Âµl O2/h/ind

%Calls subfunction *mineral_fluxes* f=1/6
f=1/6;
%L0       = (Wd0/conversion/d_V/ (1 + f * w)).^(1/ 3);  % cm, structural length
L0         = L_O2_L*del_M;
pars_tj    = [g; k; l_T; v_Hb; v_Hp]; 
temp_pars  = T_A;
% calculate mineral fluxes for an individual of structural length L0
 J_M  = mineral_fluxes(L0, f, L_m, E_m, p_M, E_G, kap, k_J, v, p_Am, p_T, E_Hp, kap_R, eta_O, n_M, n_O, pars_tj);
         
% dioxygen consumption
   TC    = tempcorr(15+273, T_ref, temp_pars);  % -, temperature correction factor
   eJO_L   = -1e6 * J_M(3,:)' .* TC* X_gas/ 24;           % µl O2/h/ind
   

 %% subfunction *mineral_fluxes* from predict_Bolinopsis_mikado by S.Augustine
     function J_M = mineral_fluxes(L,f,L_m,E_m, p_M, E_G, kap, k_J, v, p_Am, p_T, E_Hp,kap_R, eta_O, n_M, n_O, pars)
% Output: 
%  J_M: 4-n-matrix of minteral fluxes (mol of mineral compound per day)
%
%  Inputs: 
% * L0: scalar: structure of the individual
% * f, scalar, -, functional response 
% * L_m, scalar, maximum structural length (cm)
% * [E_m], scalar, maximum reserve density (J/ cm^3)
% * [p_M], scalar, volume-linked somatic maintenance costs (J/cm^3/d)
% * [E_G], scalar, cost for structure (J/cm^3)
% * kap, scalar, allocation fraction to soma (-)
% * k_J, scalar, maturity maintenance rate coefficient (1/d)
% * v, scalar, energy conductance, cm/d
% * p_Am, scalar, surface area-specific assimilation rate (J/cm^2/d)
% * p_T, scalar, surface area-linked somatic maintenance costs
% * E_H^p, scalar, maturity level at puberty (J)
% * kap_R, scalar reproduction efficiency (-)
% * eta_O, 4-4-matrix of mass to power couplers (mol/J)
% * n_M, 4-4-matrix of chemical indices, mol of element (C, H, O and N) per
% mol of compound (CO2, H2O, O2 and NH3)
% * n_O, 4-4-matrix of chemical indices, mol of element (C, H, O and N) per
% c-mol of  generalized compound  (X, V, E or P)
% * *pars_tj*, 6-1-vector of parameters for DEBtool routine get_tj pars_tj  = [g k l_T v_Hb v_Hj v_Hp]
     
     
if length(pars) == 6
 [lj lp lb info] = get_lj(pars, f);
    if info ~= 1 % numerical procedure failed
      fprintf('warning: invalid parameter value combination for get_lj \n')
    end
    L_b = lb * L_m;  L_j  = lj   * L_m;  L_p = lp * L_m; % cm, structural lengths
    ML = [ ones(length(L(L<L_b)),1) ; L(L>=L_b & L<= L_j) /L_b; ones(length(L(L>L_j)),1) * lj/ lb]; % shape correction function
 else
     [lp lb info] = get_lp (pars, f);
        if info ~= 1 % numerical procedure failed
      fprintf('warning: invalid parameter value combination for get_lp \n')
       end
     ML = 1;   L_b = lb * L_m;  L_p = lp * L_m; % cm, structural lengths
end
 
  % -------------- powers -------------------------------
  p_S  = p_M + p_T * ML./ L;                                                                % [p_S],  J/ d/ cm^3 somatic maintenance
  p_C  = f * E_m * (E_G * v * ML./ L + p_S)./ (kap * f * E_m + E_G); % [p_C], J/d/cm^3, mobilization flux
  
  p_A  =  f * p_Am .* ML .* L.^2 .* (L> L_b); % J/ d, assimilation power
  p_G = (kap * p_C - p_S) .* L.^3;                % J/d, growth power
  
  p_D1 =  p_M * L.^3 + (1- kap) * p_C .* L.^3 .* (L <=L_p);    
  p_D2 =  (k_J * E_Hp + (1-kap_R) * ((1-kap) * p_C .* L.^3 - k_J * E_Hp)) .* (L > L_p) ;   
  p_D = p_D1 + p_D2; % J/d dissipation power
  
  J_O = eta_O * [p_A p_D p_G]';               % J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M = -(n_M\n_O) * J_O;                       % J_C, J_H, J_O, J_N in rows, A, D, G in cols
   