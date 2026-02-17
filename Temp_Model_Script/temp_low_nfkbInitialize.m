function [params, species] = temp_low_nfkbInitialize()
% Initialize parameter set and species definitions for the NF-kB ODE model.
% 34C

% Core pathway parameters
params(1,1) = 0.8275; 
params(2,1) = 1.6549e-05; 
params(3,1) = 1.6549;
params(4,1) = 14.8943; 
params(5,1) = 4.1373e-07; 
params(6,1) = 4.1373e-05; 
params(6,2) = 2.4;
params(6,3) = 0.33; 
params(6,4) = 14;
params(7,1) = 0.0478;
params(8,1) = 24.8239;
params(8,2) = 1; 
params(9,1) = 0.0186; 
params(9,2) = 3.5; 
params(10,1) = 0.4965;
params(10,2) = 3.5; 
params(11,1) = 0.1303; 
params(11,2) = 0.285714; 
params(12,1) = 0.0348; 
params(12,2) = 0.285714; 
params(13,1) = 0; 
params(13,2) = 3.5; 
params(14,1) = 0.6851;
params(14,2) = 0.285714; 
params(15,1) = 0.0637; 
params(16,1) = 0.0637;
params(17,1) = 165.4945; 
params(18,1) = 165.4945; 
params(19,1) = 0.0066;
params(20,1) = 0.0066; 
params(21,1) = 157.2198; 
params(22,1) = 157.2198; 
params(23,1) = 31.4440; 
params(24,1) = 31.4440;
params(25,1) = 1.6549; 
params(26,1) = 1.6549; 
params(27,1) = 8.2747e-09; 
params(28,1) = 7.4472e-07; 
params(28,2) = 2.4; 
params(28,3) = 0.33; 
params(28,4) = 45; 
params(29,1) = 0.0042; 
params(30,1) = 24.8242; 
params(30,2) = 1;
params(31,1) = 0.0140; 
params(31,2) = 3.5;
params(32,1) = 0.4965;
params(32,2) = 3.5;
params(33,1) = 0.1303;
params(33,2) = 0.285714;
params(34,1) = 0.0348; 
params(34,2) = 0.285714; 
params(35,1) = 0; 
params(35,2) = 3.5; 
params(36,1) = 0.3426; 
params(36,2) = 0.285714; 
params(37,1) = 0.0061; 
params(38,1) = 0.0061; 
params(39,1) = 46.2806; 
params(40,1) = 46.2806; 
params(41,1) = 2.9590e-05; 
params(42,1) = 2.9590e-05; 
params(43,1) = 38.1820; 
params(44,1) = 38.1820; 
params(45,1) = 31.4440; 
params(46,1) = 31.4440; 
params(47,1) = 1.1005; 
params(48,1) = 1.1005; 
params(49,1) = 2.4824e-07; 
params(49,2) = 2.938; 
params(49,3) = 0.1775;
params(49,4) = 45; 
params(50,1) = 0; 
params(50,2) = 3.5; 
params(51,1) = 0; 
params(51,2) = 3.5; 
params(52,1) = 0.6851; 
params(52,2) = 0.285714; 
params(53,1) = 0.3426; 
params(53,2) = 0.285714; 
params(54,1) = 103.5115; 
params(55,1) = 46.2806; 
params(56,1) = 103.5107; 
params(57,1) = 46.2801; 
params(58,1) = 0.0530; 
params(59,1) = 0.0662; 
params(60,1) = 0.0530; 
params(61,1) = 0.0662; 
params(62,1) = 157.2198; 
params(63,1) = 38.1820; 
params(64,1) = 31.4440; 
params(65,1) = 31.4440; 
params(66,1) = 1.6549; 
params(67,1) = 1.1005; 

% TNF
params(68,1) = 0.0024; 
params(69,1) = 6.8051e-06; 
params(70,1) = 0.0197; 
params(71,1) = 910.2197; 
params(71,2) = 0.001; 
params(72,1) = 0.0174; 
params(72,2) = 0.001; 
params(73,1) = 0.1034; 
params(74,1) = 28.2003; 
params(75,1) = 0.0315; 
params(76,1) = 0.1034; 
params(77,1) = 1.1082; 
params(78,1) = 265.0394;
params(79,1) = 0.1034;
params(80,1) = 1.5631e3;
params(81,1) = 0.4293; 
params(82,1) = 15.5482; 
params(82,2) = 2; 
params(82,3) = 0.001116;

% IL1b
params(83,1) = 4.1374e-06; 
params(84,1) = 0.0124; 
params(85,1) = 910.2197; 
params(85,2) = 0.001; 
params(86,1) = 0.0174; 
params(86,2) = 0.001; 
params(87,1) = 0.0827; 
params(88,1) = 33.0989; 
params(89,1) = 0.0414; 
params(90,1) = 0.0827;
params(91,1) = 1.1082;
params(92,1) = 248.2417;
params(93,1) = 0.0827; 
params(94,1) = 1.4067e3;


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




