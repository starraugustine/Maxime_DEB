%% mydata_my_pet
% Sets referenced data

%%
function [data, auxData, metaData, txtData, weights] = mydata_Oikopleura_dioica  
  % created by Starrlight Augustine, Bas Kooijman, Dina Lika, Goncalo Marques and Laure Pecquerie 2015/03/31

  %% Syntax
  % [data, txt_data, metaData] = <../mydata_my_pet.m *mydata_my_pet*>
  
  %% Description
  % Sets data, pseudodata, metaData, explanatory text, weight coefficients.
  % Meant to be a template in add_my_pet
  %
  % Output
  %
  % * data: structure with data
  % * txt_data: text vector for the presentation of results
  % * metaData: structure with info about this entry
 
%% set metaData

metaData.phylum     = 'Chordata'; 
metaData.class      = 'Appendicularia'; 
metaData.order      = 'Copelata'; 
metaData.family     = 'Oikopleuridae';
metaData.species    = 'Oikopleura_dioica'; 
metaData.species_en = 'Sea squirts'; 
metaData.T_typical  = C2K(15); % K
metaData.data_0     = {'ab'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wdb'; 'Wdp'; 'Wdi'; 'Ni'};                            % tags for different types of zero-variate data
metaData.data_1     = {'t-L_fT ','t-R_fT', 'X-pX', 't-JO_f'}; % tags for different types of uni-variate data

metaData.COMPLETE   = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Starrlight Augustine','Maxime Vaugeois'};    % put names as authors as separate strings:  {'author1','author2'} , with corresponding author in first place 
metaData.date_subm = [2015 05 12];                                  % [year month day], date of entry is accepted into collection
metaData.email    = {'maxime.vaugeois@outlook.com'};               % e-mail of corresponding author
metaData.address  = {'Aix-Marseille University, 13009, France'};   % affiliation, postcode, country of the corresponding author

% uncomment and fill in the following fields when the entry is updated:
% metaData.author_mod_1  = {'author2'};                       % put names as authors as separate strings:  {'author1','author2'} , with corresponding author in first place 
% metaData.date_mod_1    = [2017 09 18];                      % [year month day], date modified entry is accepted into the collection
% metaData.email_mod_1   = {'myname@myuniv.univ'};            % e-mail of corresponding author
% metaData.address_mod_1 = {'affiliation, zipcode, country'}; % affiliation, postcode, country of the corresponding author

%% set data
% zero-variate data;
% typically depend on scaled functional response f.
% here assumed to be equal for all real data; the value of f is specified in pars_init_my_pet.

% age 0 is at onset of embryo development  
data.ab = 0.55;      units.ab = 'd';    label.ab = 'age at birth';      bibkey.ab = 'Fenaux1998';    %first house inflate ; this value is higher in Troedsson (15.7/24=0.65)  
   temp.ab = C2K(15);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
  % observed age at birth is frequently larger than ab, because of diapauzes during incubation 
data.ap = 2.7;       units.ap = 'd';    label.ap = 'age at puberty';    bibkey.ap = 'Troedsson2002'; %STA-15 fig 3(a) at 65h 
   temp.ap = C2K(15);  units.temp.ap = 'K'; label.temp.ap = 'temperature';
  % observed age at puberty is frequently larger than ap, 
  %   because allocation to reproduction starts before first eggs appear
data.am = 6.6250;    units.am = 'd';    label.am = 'life span';         bibkey.am = 'Troedsson2002'; %STA-15 tatLe 3
temp.am = C2K(15);  units.temp.am = 'K'; label.temp.am = 'temperature';
% (accounting for aging only) 

% Please specify what type of length measurement is used for your species.
% We put here snout-to-vent length, but you should change this depending on your species:
data.Lb  = 0.0115;   units.Lb  = 'cm';   label.Lb  = 'total body length at birth';  bibkey.Lb  = 'Fenaux1998';     %0.0115 is total body length but gonad is nearly null
data.Lp  = 0.0256;   units.Lp  = 'cm';   label.Lp  = 'trunk length at puberty';     bibkey.Lp  = 'Troedsson2002';  % STA-15 at 65h , tL=256.08
data.Li  = 0.0626;   units.Li  = 'cm';   label.Li  = 'ultimate trunk length';       bibkey.Li  = 'Troedsson2002';  %STA-15 at 158h Body L = 921.31 and TL/tL = 0.68 
data.Wdb = 0.015;    units.Wdb = 'mugC'; label.Wdb = 'dry weight at birth';         bibkey.Wdb = 'Troedsson2002';       %here again, this egg weight = 0.015 from Troedsson (a very lower value of 0.037 is reported in Lombard 2009)
data.Wdp = 0.36;     units.Wdp = 'mugC'; label.Wdp = 'dry weight at puberty';       bibkey.Wdp = 'Troedsson2002';  %STA-15 tatLe 2 (a), at 65h, Body weight is 0.8 µg AFDW eq to 256.08 tL ; Cw is XpXumed as 0.45 AFDW in Nakamura 1997 (so 0.36 µgC) ; Cw=10^(-6.84)*tL^2.59 in Lopez-urrutia 2003 (so 0.24µC)
data.Wdi = 6.48;     units.Wdi = 'mugC'; label.Wdi = 'ultimate dry weight';         bibkey.Wdi = 'Troedsson2002';  %-STA-15 tatLe 2 (a), at 158h, Body weight is 14.5 µg AFDW eq to 921.31 tL ; Cw is XpXumed as 0.45 AFDW in Nakamura 1997 (so 6.48µgC) ; Cw=10^(-6.84)*tL^2.59 in Lopez-urrutia 2003 (so 6.88µgC)
data.Ni  = 303;      units.Ni  = '#';    label.Ni  = 'number of eggs at death ';    bibkey.Ni  = 'Troedsson2002';  %note that a very lower value is reported in Lombard 2009 : 164
  % for an individual of ultimate length Li 
 temp.Ni = C2K(15);  units.temp.Ni = 'K'; label.temp.Ni = 'temperature';
auxData.tdeath.Ni = 6.625; units.tdeath.Ni = 'day'; label.tdeath.Ni = 'time since birth at death, Troe2002'; 

% uni-variate data

% Gonad and trunk length Troedsson2002

% at f = 0.89 (this value should be added in pars_init_my_pet as a parameter f1_tLT) 
data.tLR_S15 = [0.34      0.55      0.8      1.13     1.62    1.84      2.13      2.5     2.75      3       3.55     3.71      4.05     4.47     4.8      5.09      5.4      5.72     6.06      6.39    6.71       7    ;  % d, time since birth
                0.00003   0.00003  0.00003   0.00002  0.0004  0.00004  0.00002  0.000029  0.0011  0.0015   0.0054    0.0046    0.0073   0.0144  0.0200    0.0165    0.0226   0.0294   0.0309    0.0342  0.0253   0.0375]'; % cm, Gonad length at f and T
units.tLR_S15 = {'d', 'cm'};  label.tLR_S15 = {'time since birth', 'gonad length'};    bibkey.tLR_S15 = 'Troedsson2002';
  temp.tLR_S15 = C2K(15);  units.temp.tLR_S15 = 'K'; label.temp.tLR_S15 = 'temperature';
 
data.tL_S15 = [0.34      0.55      0.8      1.13     1.62    1.84      2.13      2.5     2.75      3       3.55     3.71      4.05     4.47     4.8      5.09      5.4      5.72     6.06      6.39    6.71       7    ;  % d, time since birth
                0.0112    0.0112   0.0131    0.0147   0.0184  0.0204    0.0230    0.0258  0.0260  0.0268   0.0339    0.0401    0.0418   0.0458  0.0534    0.0419    0.0410   0.0614   0.0548    0.0583  0.0564   0.0613  ]'; % cm, Trunk length at f and T
units.tL_S15 = {'d', 'cm'};  label.tL_S15 = {'time since birth', 'gonad length'};    bibkey.tL_S15 = 'Troedsson2002';
  temp.tL_S15 = C2K(15);  units.temp.tL_S15 = 'K'; label.temp.tL_S15 = 'temperature';
  
%at f = 0.58 (this value should be added in pars_init_my_pet as a parameter f2_tLT) 
data.tLR_L15 = [0.3324878     0.5467390     0.7905711     1.1230155     1.6180677     1.8323189     2.1204789     2.5416657     2.7412133     3.0588673     3.5467194     3.7167007     4.4708091     4.7888534     5.0917604     5.4162385     5.7136226     6.0607996     6.3785114     6.7183438     6.9621758     7.2874347     7.6272526 ; % d, time since birth
                1.955908e-05  1.905773e-05  2.958699e-05  1.115529e-05  3.685831e-05  3.622030e-05  1.680363e-05  1.896895e-05  6.642425e-04  6.923640e-04  2.901758e-03  2.779891e-03  5.670514e-03  1.015929e-02  1.086069e-02  6.933703e-03  1.585324e-02  1.792870e-02  1.668215e-02  1.598216e-02  1.745365e-02  2.158373e-02  1.548377e-02]'; % cm,Gonad length at f and T
units.tLR_L15 = {'d', 'cm'};           label.tLR_L15 = {'time since birth', 'gonad length'};   bibkey.tLR_L15 = 'Troedsson2002';
  temp.tLR_L15 = C2K(15);  units.temp.tLR_L15 = 'K'; label.temp.tLR_L15 = 'temperature';
 
data.tL_L15 = [0.3324878     0.5467390     0.7905711     1.1230155     1.6180677     1.8323189     2.1204789     2.5416657     2.7412133     3.0588673     3.5467194     3.7167007     4.4708091     4.7888534     5.0917604     5.4162385     5.7136226     6.0607996     6.3785114     6.7183438     6.9621758     7.2874347     7.6272526 ; % d, time since birth
                0.008400081   0.008302150   0.013104835   0.010969755   0.017306595   0.017273941   0.019455284   0.023732892   0.026465636   0.025130880   0.031787012   0.033553299   0.042425310   0.043934108   0.043224070   0.036604165   0.054046055   0.051356987   0.052901738   0.053443073   0.052270330   0.050060890   0.055863019 ]' ; % cm,Trunk length at f and T
units.tL_L15 = {'d', 'cm'}; label.tL_L15 = {'time since birth', 'trunk length'};   bibkey.tL_L15 = 'Troedsson2002';
  temp.tL_L15 = C2K(15);  units.temp.tL_L15 = 'K'; label.temp.tL_L15 = 'temperature';
  
  
% --------- body length at 20 deg (trunk length + gonad length)
 
%  f = 0.89 (this value should be added in pars_init_my_pet as a parameter f1_tLT)          
data.tLtot_S20  = [0.3526288   0.5901847   1.1012231   1.2667385   1.5185761   1.9215069   2.2590580   2.6038787   2.7949372   3.0816479   3.4136619   3.8290950   4.1317496   4.4250744   4.7411912   5.0714493   5.4093750   ; % d, time since birth
                0.01769144  0.01887846  0.02346265  0.02564068  0.02930905  0.03422993  0.0462147   0.0541170   0.0774069   0.0839680   0.0797031   0.0899602   0.0887187   0.0965302   0.0990892   0.1024948   0.1040988]'; % cm,Body length at f and T
units.tLtot_S20 = {'d', 'cm'}; label.tLtot_S20 = {'time since birth', 'body length'};  bibkey.tLtot_S20 = 'Troedsson2002';
  temp.tLtot_S20 = C2K(20);  units.temp.tLtot_S20 = 'K'; label.temp.tLtot_S20 = 'temperature';

%at f = 0.58 (this value should be added in pars_init_my_pet as a parameter f2_tLT)               
data.tLtot_L20  = [0.3455230  0.5974777   1.1012700   1.2667853   1.5186229   1.9289404   2.2666320   2.5971944   2.8054612   3.1076475   3.4167289   3.8409067   4.1286709   4.4306231   4.7475361   5.0857193   5.4244645   5.9212681   6.3392766; % d, time since birth
                0.0115292  0.0136109   0.0218526   0.0242454   0.0281795   0.0293765   0.0413933   0.0489085   0.0541303   0.0558680   0.0600570   0.0645388   0.0660501   0.0689147   0.0679994   0.0692506   0.0667619   0.0664631   0.0632938]'; % cm,Body length at f and T
units.tLtot_L20 = {'d', 'cm'};  label.tLtot_L20 = {'time since birth', 'body length'}; bibkey.tLtot_L20 = 'Troedsson2002';
  temp.tL_L20 = C2K(20);  units.temp.tL_L20={'K'}; label.temp.tL_L20 = {'temperature'};
  
%% Assimilation rate for different food concentration  Lombard2009a         
data.XpA     = [0.823870  1.114862  1.187610  1.933276  2.078772  2.078772  4.006592  4.443078  6.771011  8.880704  10.753962  12.609034; % XpXimilation data from Lombard 2009, Temp =15, u=gC/ind/d --> ind 784 m (tL) estim. 584 (TL) 
               0.175632  0.118343  0.129801  0.106885  0.149443  0.182179  0.146169  0.065964  0.041412  0.085606  0.080696   0.151079 ]'; % conversion factor 32e-3 J/mugC             
units.XpA = {'J/l', 'J/ind/d'}; label.XpA = {'assimilation rate', 'food concentration',};  bibkey.XpA = 'Lombard2009a';
temp.XpA = C2K(15);  units.temp.XpA = 'J/d'; label.temp.XpA = 'temperature';
auxData.Lw.XpA = 0.0584; units.Lw.XpA = 'cm'; label.Lw.XpA = 'trunk length';     


%% Lomb2005
% Respiration (O2 conso) and growth (both for gonad and trunk)  for two
% food levels

% f1 = 0.89 (this value is included in pars_init_my_pet as parameter f1_Lomb2005)

data.tJO_f1  = [1.00000   2.00000   3.00000   3.99329   4.99329;   %d, time since birth
                0.00372999	0.00977461	0.0142763	0.04443714	0.07922822]'; %micro l O2/ ind/h, O2 conso
units.tJO_f1 = {'d', 'micro l O2/ ind/h'};  label.tJO_f1 = {'time since birth', 'O2 consomation'}; bibkey.tJO_f1 = 'Lombard2005';
temp.tJO_f1 = C2K(15);  units.temp.tJO_f1 = 'K'; label.temp.tJO_f1 = 'temperature';
  
data.tL_f1  = [1.00000   2.00000   3.00000   3.99329   4.99329;   %d, time since birth
               230.0724638	266.3043478	442.0289855	668.4782609	815.2173913]'; %mum, body length
units.tL_f1 = {'d', 'mum'};  label.tL_f1 = {'time since birth', 'body length'}; bibkey.tL_f1 = 'Lombard2005';
temp.tL_f1 = C2K(15);  units.temp.tL_f1 = 'K'; label.temp.tL_f1 = 'temperature';

data.tLR_f1  = [1.00000   2.00000   3.00000   3.99329   4.99329;   %d, time since birth
                1.8115942	7.24637681	70.65217391	144.9275362	228.2608696]'; %mum, gonad length
units.tLR_f1 = {'d', 'mum'};  label.tLR_f1 = {'time since birth', 'gonad length'}; bibkey.tLR_f1 = 'Lombard2005';
temp.tLR_f1 = C2K(15);  units.temp.tLR_f1 = 'K'; label.temp.tLR_f1 = 'temperature';

% f2 = 0.58 (this value is included in pars_init_my_pet as parameter f2_Lomb2005)

data.tJO_f2  = [1.00000   2.00000   3.00000   3.99329   4.99329;   %d, time since birth
                0.00372999	0.00977461	0.01022437	0.02533699	0.04469398]'; %micro l O2/ ind/h, O2 conso
units.tJO_f2 = {'d', 'micro l O2/ ind/h'};  label.tJO_f2 = {'time since birth', 'O2 consomation'}; bibkey.tJO_f2 = 'Lombard2005';
temp.tJO_f1 = C2K(15);  units.temp.tJO_f1 = 'K'; label.temp.tJO_f1 = 'temperature';
  
data.tL_f2  = [1.00000   2.00000   3.00000   3.99329   4.99329;   %d, time since birth
               230.0724638	248.1884058	416.6666667	570.6521739	663.0434783]'; %mum, body length
units.tL_f2 = {'d', 'mum'};  label.tL_f2 = {'time since birth', 'body length'}; bibkey.tL_f2 = 'Lombard2005';
temp.tL_f2 = C2K(15);  units.temp.tL_f2 = 'K'; label.temp.tL_f2 = 'temperature';

data.tLR_f2  = [1.00000   2.00000   3.00000   3.99329   4.99329;   %d, time since birth
                1.8115942	7.24637681	56.15942029	79.71014493	119.5652174]'; %mum, gonad length
units.tLR_f2 = {'d', 'mum'};  label.tLR_f2 = {'time since birth', 'gonad length'}; bibkey.tLR_f2 = 'Lombard2005';
temp.tLR_f2 = C2K(15);  units.temp.tLR_f2 = 'K'; label.temp.tLR_f2 = 'temperature';

%% make grouped plots

set1 = {'tL_S15','tLR_S15'}; comment1 = {'Length for body, gonad at f1, 15 deg (Troes2002)'};
set2 = {'tL_L15','tLR_L15'}; comment2 = {'Length for body, gonad at f2, 15 deg (Troes2002)'};
set3 = {'tLtot_S20','tLtot_L20'}; comment3 = {'BOdy lenght 20deg, f1, f2 (Troes2002)'};
set4 = {'tJO_f1','tJO_f2'}; comment4 = {'Respiration at f1, f2 (Lomb2005)'};
set5 = {'tL_f1','tLR_f1'}; comment5 = {'Length for body, gonad (f1, Lomb2005)'};
set6 = {'tL_f2','tLR_f2'}; comment6 = {'Length for body, gonad (f2, Lomb2005)'};

metaData.grp.sets = {set1; set2; set3; set4; set5; set6};
metaData.grp.comment = {comment1;comment2;comment3; comment4; comment5; comment6};

%% set weights for all real data
weights = setweights(data, []);

%% overwriting weights (remove these remarks after editing the file)
% the weights were set automatically with the function setweigths,
% if one wants to ovewrite one of the weights it should always present an explanation example:
%
% zero-variate data:
% weights.Wdi = 100 * weights.Wdi; % Much more confidence in the ultimate dry
%                                % weights than the other data points
% uni-variate data: 
% weights.tL = 2 * weights.tL;

%% set pseudodata and respective weights
% (pseudo data are in data.psd and weights are in weights.psd)
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% overwriting pseudodata and respective weights (remove these remarks after editing the file)
% the pseudodata and respective weights were set automatically with the function setpseudodata
% if one wants to ovewrite one of the values it should always present an explanation
% example:
% data.psd.p_M = 1000;                    % my_pet belongs to a group with high somatic maint 
% weight.psd.kap = 10 * weight.psd.kap;   % I need to give this pseudo data a higher weight

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
if exist('comment','var')
  txtData.comment = comment;
end

%% References
  %
  bibkey = 'Troedsson2002'; type = 'Article'; bib = [ ... 
  'author = {Troedsson, C. and Bouquet, J. M. and Aksnes, D. L. and Thompson, E. M.}, ' ... 
  'year = {2002}, ' ...
  'title = {Resource allocation between somatic growth and reproductive output in the pelagic chordate Oikopleura dioica allows opportunistic response to nutritional variation} ' ...
  'journal = {Marine Ecology Progress Series}, ' ...
  'volume = {243}, ' ...
  'pages = {83-91}'];
  metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

  %
  bibkey = 'Lombard2005'; type = 'Article'; bib = [ ... 
  'author = {TLombard, F. and Sciandra, A. and Gorsky, G.}, ' ... 
  'year = {2005}, ' ...
  'title = {Influence of body mXpX, food concentration, temperature and filtering activity on the oxygen uptake of the appendicularian Oikopleura dioica} ' ...
  'journal = {Marine Ecology Progress Series}, ' ...
  'volume = {301}, ' ...
  'pages = {149-158}'];
  metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

  %
  bibkey = 'Lombard2009a'; type = 'Article'; bib = [ ... 
  'author = {TLombard, F. and Sciandra, A. and Gorsky, G.}, ' ... 
  'year = {2009}, ' ...
  'title = {{Appendicularian ecophysiology II: Modeling nutrition, metabolism, growth and reproduction of the appendicularian Oikopleura dioica} ' ...
  'journal = {Journal of Marine Systems}, ' ...
  'volume = {78}, ' ...
  'pages = {606-616}'];
 metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

  %
  bibkey = 'Lombard2009b'; type = 'Article'; bib = [ ... 
  'author = {TLombard, F. and Sciandra, A. and Gorsky, G.}, ' ... 
  'year = {2009}, ' ...
  'title = {{Appendicularian ecophysiology II: Modeling nutrition, metabolism, growth and reproduction of the appendicularian Oikopleura dioica} ' ...
  'journal = {Journal of Marine Systems}, ' ...
  'volume = {78}, ' ...
  'pages = {617-629}'];
  metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

  %
  bibkey = 'Fenaux1998'; type = 'Article'; bib = [ ... 
  'author = {Fenaux, R.}, ' ... 
  'year = {1998}, ' ...
  'booktitle={The Biology of Pelagic Tunicates},' ...
  'title = {Anatomy and functional morphology of the Appendicularia} ' ...
  'pages={25-34},' ...
  'editor={Bone, Q.},' ...
  'putLisher={Oxford University Press: Oxford}'];
  metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
  %

%% Facts
% * Standard model with egg (not foetal) development and no acceleration
  
%% Discussion points
pt1 = 'Author_mod_1: ';
pt2 = 'Author_mod_1: ';     
metaData.discussion = {pt1; pt2}; 

