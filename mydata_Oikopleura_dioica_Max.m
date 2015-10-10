%% mydata_Oikopleura_dioica
close all
clear all
clc

addpath(genpath('../DEBtool_O/'))

COMPLETE = 2.4;   
FIT = 8.1897;         
                      
  global n_O n_M dwm ass_par CONC X_gas L_O2 L_O2_L % pass d_O, w_O, mu_O directly to predict_Oikopleura_dioica

%% set data

% real data
ab = 0.6;           %  1 d, age at birth at f (age 0 is at onset of embryo development) from Fenaux,1998 (0.55); Troedsson (15.7/24=0.65).
  T_ab = 273 + 15;       % K, temperature for ab. from Fenaux,1998 and Troedsson2002.
ap = 2.75;          %  2 d, age at puberty at f from Troedsson2002-STA-15 fig 3(a) at 66h
  T_ap = 273 + 15;       % K, temperature for ap. from Troedsson2002 
Lb = 0.0115;        %  3 cm, physical length at birth at f. from Fenaux,1998 - 0.0115 is total body length but gonad is nearly null
Lp = 0.0256;        %  4 cm, physical length at puberty at f: in Troedsson2002-STA-15 at 66h , BL=256.08
L8 = 0.0626;        %  5 cm, physical length at death (6.6 d) at f. from Troedsson2002-STA-15 at 158h BL = 921.31 and TL/BL = 0.68 
Wb = 0.037;         %  6 mugC, Carbon weight at birth at f. egg weight = 0.015 from Troedsson and 0.037 from Lombard
Wp = 0.36;          %  7 mugC, Carbon weight at puberty at f: in Troedsson2002-STA-15 table 2 (a), at 66h, Body weight is 0.8 µg AFDW eq to 256.08 BL ; Cw is assumed as 0.45 AFDW in Nakamura 1997 (so 0.36 µgC) ; Cw=10^(-6.84)*BL^2.59 in Lopez-urrutia 2003 (so 0.24µC)
W8 = 6.48;          %  8 mugC, Carbon weight at death (6.6 d) at f: Troedsson2002-STA-15 table 2 (a), at 158h, Body weight is 14.5 µg AFDW eq to 921.31 BL ; Cw is assumed as 0.45 AFDW in Nakamura 1997 (so 6.48µgC) ; Cw=10^(-6.84)*BL^2.59 in Lopez-urrutia 2003 (so 6.88µgC)
NR = 164;           %  9 #, number of eggs at death (6.6 d) at f from Troedsson2002 =303 (note that a very lower value is reported in Lombard 2009 : 164)
  T_R = 273 + 15;        % K, temperature for Ri
tm = 6.6250;        % 10 d, life span at f from Troedsson2002-STA-15 table 3, but this value does not relate to ageing
  T_tm = 273 + 15;       % K, temperature for tm

% pseudo-data from pars_my_pet at T_ref; don't change these data
v = 0.02;           % 11 cm/d, energy conductance
kap = 0.8;          % 12 -, alloaction fraction to soma = growth + somatic maintenance
kap_R = 0.95;       % 13 -, reproduction efficiency
p_M = 18;           % 14 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T =  0;           % 15 J/d.cm^2, {p_T} surface-specific som maintenance
k_J = 0.001;        % 16 1/d, < k_M = p_M/E_G, maturity maint rate coefficient
kap_G = .8;         % 17 -, growth efficiency

% pack data
data = [ab; ap; Lb; Lp; L8; Wb; Wp; W8; NR; tm; % 01:10 real data
        v; kap; kap_R; p_M; p_T; k_J; kap_G];   % 11:16 pseudo data

% weight coefficients for WLS
%Data = [data(:,[1 1]), min(100,1 ./ max(1e-6, data) .^ 2)]; % nmregr, nrregr (WLS criterion)
%Data(2,3)=Data(2,3)*15;
%Data(8,3)=Data(8,3)*200;
%Data(9,3)=Data(9,3)*1000;
%Data(10,3)=Data(10,3)*10;
%%Data(1:2,3) = 100./ data(1:2) .^ 2;   
%%Data(3:5,3) = 1./ data(3:5) .^ 2; 
%%Data(6:8,3) =  1./ data(6:8) .^ 2;                              % give age at birt data weight
%%Data(8:10,3) =  100./ data(8:10) .^ 2;
%Data(11:16,3) = 0;  
%Data(17,3) = 100./ data(17) .^ 2                               % more weight to kap_G
%%
Data = [data(:,[1 1]), min(100,1 ./ max(1e-6, data) .^ 2)]; % nmregr, nrregr (WLS criterion)
Data(1:9,3)      = 10 * Data(1:9,3);       % give real data more weight
Data([6 8],3)      = 10 * Data([6 8],3);      % give weight data more weight
Data([3 5],3)      = 100 * Data([3 5],3);      % give length data more weight
Data([1 2],3)  = 3e4 * Data([1 2],3);   % give age data more weight
Data([11 16],3)  = 1e-3 * Data([11 16],3);   % give pseudo-data less weight
Data(17,3)        = 1e5 * Data(17,3);        % more weight to kap_G
Data(5,3)        = 1e4 * Data(5,3);        % more weight to physical length at death
Data(7,3)        = 5e3 * Data(7,3);        % more weight to Carbon weight at death 
Data(8,3)        = 5e3 * Data(8,3);        % more weight to Carbon weight at death 
Data(9,3)        = 1e-2 * Data(9,3);        % more weight to NR
Data(2,3)        = 1e1 * Data(2,3);        % more weight to age at puberty
Data(10,3)        = 1e2 * Data(2,3);        % more weight to age at puberty
% insert temperature data for rates and times in the first column
Data([1 2 9 10], 1) = [T_ab; T_ap; T_R; T_tm];

txt_data = {... % for presentation of predictions
    '1 ab, d, age at birth ';
    '2 ap, d, age at puberty ';
    '3 Lb, cm, physical length at birth ';
    '4 Lp, cm, physical length at puberty ';
    '5 L8, cm, physical length at death';
    '6 Wb, mug, Carbon weight at birth ';
    '7 Wp, mug, Carbon weight at puberty ';
    '8 W8, mug, Carbon weight at death ';
    '9 NR, #, number of eggs at death ';
   '10 am, d, life span ';
   '11 v, cm/d, energy conductance ';
   '12 kap, -, allocation fraction to soma ';
   '13 kap_R, -, reproduction efficiency ';
   '14 [p_M], J/d.cm^3, vol-spec som maint ';
   '15 {p_T}, J/d.cm^2, sur-spec som maint ';
   '16 k_J, 1/d, maturity maint rate coefficient ';
   '17 kap_G, -, growth efficiency'};


%  This info is copied directly in predict_my_pat

%tgL_S15 --> Troedsson 2002, Standard (f=1), T=15 + 273, trunk & gonad
 tgL_S15 = [0.34      0.55      0.8      1.13     1.62    1.84      2.13      2.5     2.75      3       3.55     3.71      4.05     4.47     4.8      5.09      5.4      5.72     6.06      6.39    6.71       7    ; % d, time since birth
      0.0112    0.0112   0.0131    0.0147   0.0184  0.0204    0.0230    0.0258  0.0260  0.0268   0.0339    0.0401    0.0418   0.0458  0.0534    0.0419    0.0410   0.0614   0.0548    0.0583  0.0564   0.0613  ; % cm,Trunk length at f and T
      0.00003   0.00003  0.00003   0.00002  0.0004  0.00004  0.00002  0.000029  0.0011  0.0015   0.0054    0.0046    0.0073   0.0144  0.0200    0.0165    0.0226   0.0294   0.0309    0.0342  0.0253   0.0375]'; % cm,Gonad length at f and T
      
      tgL_S15 = [tgL_S15, 10 ./ tgL_S15(:,2:3).^2]; % append weight coefficients for WLS criterion

tL_S15  = tgL_S15(:,[1 2 4]);
gL_S15  = tgL_S15(:,[1 3 5]); 

%tgL_L15 --> Troedsson 2002, Limited (f=1/2), T=15 + 273, trunk & gonad
 tgL_L15 = [0.3324878 0.5467390 0.7905711 1.1230155 1.6180677 1.8323189 2.1204789 2.5416657 2.7412133 3.0588673 3.5467194 3.7167007 4.4708091 4.7888534 5.0917604 5.4162385 5.7136226 6.0607996 6.3785114 6.7183438 6.9621758 7.2874347 7.6272526 ; % d, time since birth
      0.008400081 0.008302150 0.013104835 0.010969755 0.017306595 0.017273941 0.019455284 0.023732892 0.026465636 0.025130880 0.031787012 0.033553299 0.042425310 0.043934108 0.043224070 0.036604165 0.054046055 0.051356987 0.052901738 0.053443073 0.052270330 0.050060890 0.055863019  ; % cm,Trunk length at f and T
      1.955908e-05 1.905773e-05 2.958699e-05 1.115529e-05 3.685831e-05 3.622030e-05 1.680363e-05 1.896895e-05 6.642425e-04 6.923640e-04 2.901758e-03 2.779891e-03 5.670514e-03 1.015929e-02 1.086069e-02 6.933703e-03 1.585324e-02 1.792870e-02 1.668215e-02 1.598216e-02 1.745365e-02 2.158373e-02 1.548377e-02]'; % cm,Gonad length at f and T
      
      tgL_L15 = [tgL_L15, 10 ./ tgL_L15(:,2:3).^2]; % append weight coefficients for WLS criterion

tL_L15  = tgL_L15(:,[1 2 4]);
gL_L15  = tgL_L15(:,[1 3 5]); 

%bL_S20 --> Troedsson 2002, Standard (f=1), T=20 + 273, total body
 bL_S20 = [0.3526288 0.5901847 1.1012231 1.2667385 1.5185761 1.9215069 2.2590580 2.6038787 2.7949372 3.0816479 3.4136619 3.8290950 4.1317496 4.4250744 4.7411912 5.0714493 5.4093750   ; % d, time since birth
      0.01769144 0.01887846 0.02346265 0.02564068 0.02930905 0.03422993 0.04621473 0.05411701 0.07740668 0.08396804 0.07970307 0.08996025 0.08871870 0.09653017 0.09908925 0.10249477 0.10409879]'; % cm,Body length at f and T
      
      bL_S20 = [bL_S20, 10 ./ bL_S20(:,2).^2];      % append weight coefficients for WLS criterion

%bL_L20 --> Troedsson 2002, Limited (f=1/2), T=20 + 273, total body
 bL_L20 = [0.3455230 0.5974777 1.1012700 1.2667853 1.5186229 1.9289404 2.2666320 2.5971944 2.8054612 3.1076475 3.4167289 3.8409067 4.1286709 4.4306231 4.7475361 5.0857193 5.4244645 5.9212681 6.3392766; % d, time since birth
      0.01152923 0.01361088 0.02185265 0.02424538 0.02817951 0.02937648 0.04139326 0.04890853 0.05413031 0.05586802 0.06005703 0.06453878 0.06605009 0.06891474 0.06799941 0.06925064 0.06676193 0.06646310 0.06329387]'; % cm,Body length at f and T
      
      bL_L20 = [bL_L20, 10 ./ bL_L20(:,2).^2];      % append weight coefficients for WLS criterion

%Assimilation data from Lombard 2009 A      
Ass = [5.488496 3.698216 4.056272 3.340160 4.670081 5.693099 4.567781 2.061388 1.294123 2.675197 2.521744 4.721234; % assimilation data from Lombard 2009, Temp =15, u=gC/ind/d --> ind 784 m (BL) estim. 584 (TL) 
        25.74593  34.83943  37.11280  60.41487  64.96162  64.96162 125.20600 138.84620 211.59410 277.52200 336.06130 394.03230]'; %conetration gC/l
CONC= Ass(:,2);
Ass = [Ass, 10 ./ Ass(:,1).^2];

%Respiration rate + Growth vs time from Lombard 2005
%          time, O2 conso (micro litre O2/ ind/h), Trunk length (cm), Gonad length (cm)
TAB    = [1,	4.15094340e-003, 2.30563924e-002, 0;
           2,	9.86937591e-003, 2.67248676e-002, 4.85962699e-04;
           3,	1.43686502e-002, 4.42419014e-002, 6.77834458e-03;
           3.99328859,	4.42962761e-002, 6.71525213e-002,1.43828420e-002;
           4.99328859, 7.91003400e-002, 8.18998758e-002 , 2.27162023e-002];
O2_H    = TAB(:,2);
L_O2    = TAB(:,3);
           
O2_H  = [O2_H , 10 ./ O2_H(:,1).^2];

TAB2    = [1,	4.15094340e-003, 2.30563924e-002, 0;
           2,	9.86937591e-003, 2.45382787e-002, 4.85962699e-04;
           3,	1.00147086e-002, 4.16179947e-002, 5.46655338e-03;
           3.99328859,	2.51380758e-002, 5.75315300e-002, 7.96884780e-002;
           4.99328859, 4.44414140e-002, 6.65937533e-002 , 1.20748029e-002];
O2_L    = TAB2(:,2);
L_O2_L    = TAB2(:,3);
           
O2_L  = [O2_L , 10 ./ O2_L(:,1).^2];

%% conversion coefficients (selected copy-paste from pars_my_pet)

% chemical indices
%       X     V     E     P
n_O = [1.00, 1.00, 1.00, 1.00;          % C/C, equals 1 by definition
        1.80, 1.80, 1.80, 1.80;         % H/C  these values show that we consider dry-mass
        0.50, 0.50, 0.50, 0.50;         % O/C
        0.15, 0.15, 0.15, 0.15];        % N/C  
            
% minerals
%   rows: elements carbon, hydrogen, oxygen, nitrogen 
%   columns: carbon dioxide (C), water (H), dioxygen (O), ammonia (N)
%     CO2 H2O O2 NH3
n_M = [1,  0, 0,  0;  % C
            0,  2, 0,  3;  % H
            2,  1, 2,  0;  % O
            0,  0, 0,  1]; % N
            
 % specific densities
%       X     V     E     P
%d_O = [0.3;  0.045;  0.045;  0.3];     % g/cm^3, specific densities for organics
d_O  = [0.03284;  0.076; 0.0695;  0.3]; % new values
%       X     V     E     P
d_OC = [0.1645;  0.038; 0.0349;  0.3];  % g C/cm^3, specific densities for organics
% chemical potentials
%       X     V     E     P
mu_O = [525; 500;  550;  480] * 1000;   % J/mol C, chemical potentials for organics X-->Brown1991+Troedsson2005

% molecular weights
w_O = n_O' * [12; 1; 16; 14];           % g/mol, mol-weights for organics

% Assimilation parameters
ass_par = [32e-3  15 (0.0784-0.02)];    % conversion factor u=J/gC , Temperature in Celsuis, Trunk length of individual during experiment in cm

% pack coefficients
dwm = [d_O, w_O, mu_O, d_OC];           % g/cm^3, g/mol, kJ/mol spec density, mol weight, chem pot

%% parameters: initial values at T_ref
T_ref  =2.9300e+002  ;  %  1 K, temp for which rate pars are given; don't change this vulue
T_A    =1.1887e+004  ;  %  2 K, Arrhenius temp 
f      =1.0000e+000  ;  %  3 -, scaled functional response
z      =2.3255e-002  ;  %  4 -, zoom factor; for z = 1: L_m = 1 cm
del_M  =3.1716e-001  ;  %  5 -, shape coefficient
F_m    =2.4000e+003  ;  %  6 l/d.cm^2, {F_m} max spec searching rate
kap_X  =8.0000e-001  ;  %  7 -, digestion efficiency of food to reserve
v      =2.7868e-002  ;  %  8 cm/d, energy conductance
kap    =5.3998e-001  ;  %  9 -, alloaction fraction to soma = growth + somatic maintenance
kap_R  =9.5000e-001  ;  % 10 -, reproduction efficiency
p_M    =1.3238e+004  ;  % 11 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T    =0.0000e+000  ;  % 12 J/d.cm^2, {p_T} surface-specific som maintenance
k_J    =1.0000e-003  ;  % 13 1/d, maturity maint rate coefficient
E_G    =2.8911e+003  ;  % 14 J/cm^3, [E_G], spec cost for structure
E_Hb   =4.2422e-005  ;  % 15 J, E_H^b maturity threshold at birth
E_Hp   =5.4941e-003  ;  % 16 J, E_H^p maturity threshold at puberty
h_a    =1.7000e-018  ;  % 17 1/d^2, Weibull aging acceleration
s_G    =1.0000e+001  ;  % 18 -, Gompertz stress coefficient
v_R    =5.7567e-002  ;  % 19 cm/d, contribution of reprod buffer to total trunc length
kap_X_P= 0.1;           %  20 -, defecation efficiency
%
T_ref   =2.9300e+002    ;  %  1 K, temp for which rate pars are given; don't change this vulue
T_A     =8000           ;  %  2 K, Arrhenius temp 
f       =1.0000e+000    ;  %  3 -, scaled functional response
z       =1.9000e-002 *1.18   ;  %  4 -, zoom factor; for z = 1: L_m = 1 cm
del_M   =2.8716e-001 *.9   ;  %  5 -, shape coefficient
F_m     =2.4000e+003    ;  %  6 l/d.cm^2, {F_m} max spec searching rate
kap_X   =8.0000e-001    ;  %  7 -, digestion efficiency of food to reserve
v       =2.7868e-002*.8  ;  %  8 cm/d, energy conductance
kap     =1.8998e-001    ;  %  9 -, alloaction fraction to soma = growth + somatic maintenance
kap_R   =9.5000e-001    ;  % 10 -, reproduction efficiency
p_M     =1.3238e+004    ;  % 11 J/d.cm^3, [p_M] vol-specific somatic maintenance
p_T     =0.0000e+000    ;  % 12 J/d.cm^2, {p_T} surface-specific som maintenance
k_J     =1.0000e-003    ;  % 13 1/d, maturity maint rate coefficient
E_G     =2.8911e+003    ;  % 14 J/cm^3, [E_G], spec cost for structure
E_Hb    =4.2422e-005*12  ;  % 15 J, E_H^b maturity threshold at birth
E_Hp    =5.4941e-003*2  ;  % 16 J, E_H^p maturity threshold at puberty
h_a     =1.7000e-018 *10000000  ;  % 17 1/d^2, Weibull aging acceleration
s_G     =6.5356   ;  % 18 -, Gompertz stress coefficient
v_R     =3.1e-002       ;  % 19 cm/d, contribution of reprod buffer to total trunc length
kap_X_P =0.25           ;  % 20 -, defecation efficiency
%

%T_ref=   2.9300e+02   
%T_A=   8.0000e+03 
%f=   1.0000e+00   
%z=   5.7770e-02   
%del_M=   1.4730e-01   
%F_m =   2.8506e+03   
%kap_X=   8.0000e-01  
%v=   2.1705e-02   
%kap=   1.8990e-03 
%kap_R =   9.5000e-01   
%p_M =   3.0508e+02   
%p_T =   0.0000e+00  
%k_J =  1.0000e-03  
%E_G =  1.8849e+03   
%E_Hb =   3.9943e-04*6
%E_Hp =  6.5074e-02   
%h_a  = 1.7175e-11   
%s_G  = 2.0179e+01   
%v_R  = 3.1e-002   
%kap_X_P=   5.6379e-02   
%   
   % molar volume of gas at 1 bar and 0 C is 22.414 L/mol
T        = 273 + 15;     
X_gas    = T*22.414/273;   % L/mol

%pack parameters and fix T_ref and f and possibly other at well
%   in second column: 0 = fix; 1 = release
pars = [T_ref 0; T_A 0; f    0; z     1; del_M 1; F_m  0;  
        kap_X 0; v   1; kap  1; kap_R 0; p_M   1; p_T  0;  
        k_J   0; E_G 1; E_Hb 1; E_Hp  1; h_a   1; s_G  1;
        v_R   0; kap_X_P 1]; 

        
txt_pars = { ...    % for presentation of parameter estimates
  'T_ref, K'; 'T_A, K';        'f, -'; 
  'z, -';     'del_M, -';      '{F_m}, l/d.cm^2'; 
  'kap_X, -'; 'v, cm/d';       'kap, -'; 
  'kap_R, -'; '[p_M], J/d.cm^3'; '{p_T}, J/d.cm^2'; 
  'k_J, 1/d'; '[E_G], J/cm^3'; 'E_Hb, J'; 'E_Hp, J';  
  'h_a, 1/d^2'; 's_G, -'; 'v_R, cm/d'};

%% estimate parameters

nmregr_options('default');                     % set options for parameter estimation
nmregr_options('max_step_number',5e4);         % set options for parameter estimation 3e2
nmregr_options('max_fun_evals',2e4);           % set options for parameter estimation


pars = nmregr('predict_Oikopleura_dioica',pars, Data, tL_S15, gL_S15,tL_L15, gL_L15,bL_S20,bL_L20,Ass,O2_H,O2_L)
sd = 0 * pars(:,1);                            % initiate standard deviations

[cov cor sd ss] = pregr('predict_Oikopleura_dioica', pars, Data, tL_S15, gL_S15,tL_L15, gL_L15,bL_S20,bL_L20,Ass,O2_H,O2_L) % get standard deviation for WLS
% get FIT
Data(:,3) = 0; Data(1:10,3) = 1;               % give unit weight to real data, zero to pseudo-data
[MRE RE] = mre('predict_Oikopleura_dioica',...
 pars, Data, tL_S15, gL_S15,tL_L15,...
 gL_L15,bL_S20,bL_L20,Ass,O2_H,O2_L);
FIT = 10 * (1 - MRE)                           % get mark for goodness of fit

%% get expectations

t    = linspace(0,7.5,100)';                   % times for plotting length data
t_20 = linspace(0,6.5,100)';                   % times for plotting length data
C    = linspace(25,400,100)'; 

  [Edata, EL_S15, GL_S15,...
  EL_L15, GL_L15,...
  BL_S20,...
  BL_L20,...
  ASS,eJO,eJO_L] = predict_Oikopleura_dioica(pars(:,1), Data, t, t, ...
  t, t,...
  t_20,...
  t_20,...
  C,t,t); % notice use of first column of pars only

%% present results

 printpar(txt_data, data, Edata, 'data and expectations'); % for zero-variate data
 fprintf('\n')                                           
 printpar(txt_pars, pars, sd);
 close all
 nbfig = 7;
 Large = 0.5;
 Long  = 0.5;
 
 set(0,'defaultaxesfontsize',20)
 set (0,'defaultlinelinewidth', 2)
 set(0,'DefaultLineMarkerSize',12);

 figure(1,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig)*((1/2)/nbfig) .4-(nbfig-nbfig)*((.4)/nbfig) Large Long])  % one figure to show results of uni-variate data
 plot(tL_S15(:,1), tL_S15(:,2), '.r', t, EL_S15, 'r')
 axis([0 7 0 0.08])
 hold on
 plot(gL_S15(:,1), gL_S15(:,2), '.g', t, GL_S15, 'g')
 xlabel('time since birth, d')
 ylabel('Length, cm')
 title('Standard 15C - Trunk & Gonad Length')
 print('plots/figure1.eps','-deps')

 figure(2,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig+1)*((1/2)/nbfig) .4-(nbfig-nbfig+1)*((.4)/nbfig) Large Long])  % one figure to show results of uni-variate data
 plot(tL_L15(:,1), tL_L15(:,2), '.r', t, EL_L15, 'r')
 axis([0 7 0 0.08])
 hold on
 plot(gL_L15(:,1), gL_L15(:,2), '.g', t, GL_L15, 'g')
 xlabel('time since birth, d')
 ylabel('Length, cm')
 title('Limited 15C - Trunk & Gonad Length')
  print('plots/figure2.eps','-deps')
 
 figure(3,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig+2)*((1/2)/nbfig) .4-(nbfig-nbfig+2)*((.4)/nbfig) Large Long])  % one figure to show results of uni-variate data
 plot(bL_S20(:,1), bL_S20(:,2), '.', t, BL_S20,'b')
 axis([0 7 0 0.11])
 xlabel('time since birth, d')
 ylabel('Length, cm')
 title('Standard 20C - Total Body Length')
   print('plots/figure3.eps','-deps')
 
 figure(4,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig+3)*((1/2)/nbfig) .4-(nbfig-nbfig+3)*((.4)/nbfig) Large Long])
 title('Limited 20')
 plot(bL_L20(:,1), bL_L20(:,2), '.', t, BL_L20,'b')
 axis([0 7 0 0.11])
 xlabel('time since birth, d')
 ylabel('Length, cm')
 title('Limited 20C - Total Body Length')
   print('plots/figure4.eps','-deps')
 
 figure(5,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig+4)*((1/2)/nbfig) .4-(nbfig-nbfig+4)*((.4)/nbfig) Large Long])
 plot(Ass(:,2).*ass_par(1),Ass(:,1).*ass_par(1), '.',Ass(:,2).*ass_par(1),ASS,'b')
 xlabel('time since birth, d')
 ylabel('mu g C/d')
 title('Assimilation rate')
   print('plots/figure5.eps','-deps')
 
 
  figure(6,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig+5)*((1/2)/nbfig) .4-(nbfig-nbfig+5)*((.4)/nbfig) Large Long])
 plot(TAB(:,1),TAB(:,2), '.',TAB(:,1),eJO,'b')
 xlabel('time since birth, d')
 ylabel('mu l O2 /h')
  title('Respiration rate at f=1')
    print('plots/figure6.eps','-deps')

 figure(7,'position',get(0,'screensize')([3,4,3,4]).*[(nbfig-nbfig+6)*((1/2)/nbfig) .4-(nbfig-nbfig+6)*((.4)/nbfig) Large Long])
 plot(TAB2(:,1),TAB2(:,2), '.',TAB2(:,1),eJO_L,'b')
 xlabel('time since birth, d')
 ylabel('mu l O2 /h')
  title('Respiration rate at f=1/6')
    print('plots/figure7.eps','-deps')