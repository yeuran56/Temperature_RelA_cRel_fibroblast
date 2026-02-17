function [params, species] = temp_exlow_nfkbInitialize()
% Initialize parameter set and species definitions for the NF-kB ODE model.
% 24C


% Core pathway parameters
params(1,1) = 0.4281; 
params(2,1) = 8.5632e-06; 
params(3,1) = 0.8563; 
params(4,1) = 7.7069; 
params(5,1) = 2.1408e-07; 
params(6,1) = 2.1408e-05; 
params(6,2) = 2.4; 
params(6,3) = 0.33; 
params(6,4) = 14;
params(7,1) = 0.0247; 
params(8,1) = 12.8448; 
params(8,2) = 1;
params(9,1) = 0.0096; 
params(9,2) = 3.5; 
params(10,1) = 0.2569; 
params(10,2) = 3.5; 
params(11,1) = 0.0674; 
params(11,2) = 0.285714; 
params(12,1) = 0.0180; 
params(12,2) = 0.285714; 
params(13,1) = 0; 
params(13,2) = 3.5; 
params(14,1) = 0.3545; 
params(14,2) = 0.285714; 
params(15,1) = 0.0330; 
params(16,1) = 0.0330; 
params(17,1) = 85.6317;
params(18,1) = 85.6317; 
params(19,1) = 0.0034; 
params(20,1) = 0.0034;
params(21,1) = 81.3502; 
params(22,1) = 81.3502;
params(23,1) = 16.2700; 
params(24,1) = 16.2700;
params(25,1) = 0.8563; 
params(26,1) = 0.8563; 
params(27,1) = 4.2815e-09; 
params(28,1) = 3.8534e-07; 
params(28,2) = 2.4; 
params(28,3) = 0.33; 
params(28,4) = 45; 
params(29,1) = 0.0022; 
params(30,1) = 12.8448; 
params(30,2) = 1; 
params(31,1) = 0.0072; 
params(31,2) = 3.5; 
params(32,1) = 0.2569; 
params(32,2) = 3.5;
params(33,1) = 0.0674; 
params(33,2) = 0.285714; 
params(34,1) = 0.0180; 
params(34,2) = 0.285714; 
params(35,1) = 0; 
params(35,2) = 3.5; 
params(36,1) = 0.1773; 
params(36,2) = 0.285714; 
params(37,1) = 0.0032; 
params(38,1) = 0.0032; 
params(39,1) = 23.9470;
params(40,1) = 23.9470;
params(41,1) = 1.5311e-05; 
params(42,1) = 1.5311e-05; 
params(43,1) = 19.7565; 
params(44,1) = 19.7565; 
params(45,1) = 16.2700; 
params(46,1) = 16.2700; 
params(47,1) = 0.5695; 
params(48,1) = 0.5695; 
params(49,1) = 1.2845e-07; 
params(49,2) = 2.938;
params(49,3) = 0.1775; 
params(49,4) = 45; 
params(50,1) = 0; 
params(50,2) = 3.5; 
params(51,1) = 0; 
params(51,2) = 3.5; 
params(52,1) = 0.3545;
params(52,2) = 0.285714; 
params(53,1) = 0.1773; 
params(53,2) = 0.285714; 
params(54,1) = 53.5599;
params(55,1) = 23.9470; 
params(56,1) = 53.5599; 
params(57,1) = 23.9470; 
params(58,1) = 0.0274; 
params(59,1) = 0.0343; 
params(60,1) = 0.0274; 
params(61,1) = 0.0343; 
params(62,1) = 81.3502; 
params(63,1) = 19.7565; 
params(64,1) = 16.2700; 
params(65,1) = 16.2700; 
params(66,1) = 0.8563;
params(67,1) = 0.5695;

% TNF
params(68,1) = 0.0012;
params(69,1) = 3.5212e-06; 
params(70,1) = 0.0102; 
params(71,1) = 470.9746;
params(71,2) = 0.001; 
params(72,1) = 0.0090; 
params(72,2) = 0.001; 
params(73,1) = 0.0535; 
params(74,1) = 14.5916; 
params(75,1) = 0.0163; 
params(76,1) = 0.0535; 
params(77,1) = 0.5734; 
params(78,1) = 137.1392; 
params(79,1) = 0.0535; 
params(80,1) = 808.7918; 
params(81,1) = 0.2221; 
params(82,1) = 8.0451; 
params(82,2) = 2; 
params(82,3) = 0.001116; 

% IL1b
params(83,1) = 2.1408e-06;
params(84,1) = 0.0064; 
params(85,1) = 470.9746; 
params(85,2) = 0.001;
params(86,1) = 0.0090; 
params(86,2) = 0.001; 
params(87,1) = 0.0428;
params(88,1) = 17.1263;
params(89,1) = 0.0214; 
params(90,1) = 0.0428;
params(91,1) = 0.5734; 
params(92,1) = 128.4476; 
params(93,1) = 0.0428; 
params(94,1) = 727.8698; 



% SPECIES
species.NAMES = {...
'stim' ...1
'IkBa' ...2
'IkBan' ...3
'IkBe' ...4
'IkBen' ...5
'IkBaRelA' ...6
'IkBaRelAn' ...7
'IkBacRel' ...8
'IkBacReln' ...9
'IkBeRelA' ...10
'IkBeRelAn' ...11
'IkBecRel' ...12
'IkBecReln' ...13
'IkBat' ...14
'IkBet' ...15
'IKKIkBaRelA' ...16
'IKKIkBacRel' ...17
'IKKIkBa' ...18
'IKKIkBeRelA' ...19
'IKKIkBecRel' ...20
'IKKIkBe' ...21
'cRel' ...22
'cReln' ...23
'RelA' ...24
'RelAn' ...25
'IKK_off' ...26
'IKK' ...27
'IKK_i' ...28
'TNF' ...29
'TNFR' ...30
'TNFR_TNF' ...31
'TTR' ...32
'C1_off' ...33
'C1' ...34
'TAK1_off' ...35
'TAK1' ...36
'IL1' ...37
'IL1R' ...38
'IL1R_IL1' ...39
'MIT' ...40
'C2_off' ...41
'C2' ...42

};

