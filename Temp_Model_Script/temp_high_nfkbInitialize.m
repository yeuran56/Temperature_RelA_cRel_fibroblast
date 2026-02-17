function [params, species] = temp_high_nfkbInitialize()
% Initialize parameter set and species definitions for the NF-kB ODE model.
% 40C


% Core pathway parameters
params(1,1) = 1.2041; 
params(2,1) = 2.4082e-05; 
params(3,1) = 2.4082; 
params(4,1) = 21.6742; 
params(5,1) = 6.0206e-07;
params(6,1) = 6.0206e-05; 
params(6,2) = 2.4; 
params(6,3) = 0.33; 
params(6,4) = 14; 
params(7,1) = 0.0696; 
params(8,1) = 36.1237;
params(8,2) = 1; 
params(9,1) = 0.0271; 
params(9,2) = 3.5; 
params(10,1) = 0.7225; 
params(10,2) = 3.5;
params(11,1) = 0.1896; 
params(11,2) = 0.285714; 
params(12,1) = 0.0506; 
params(12,2) = 0.285714; 
params(13,1) = 0; 
params(13,2) = 3.5; 
params(14,1) = 0.9970; 
params(14,2) = 0.285714; 
params(15,1) = 0.0927;
params(16,1) = 0.0927; 
params(17,1) = 240.8245; 
params(18,1) = 240.8245; 
params(19,1) = 0.0096;
params(20,1) = 0.0096; 
params(21,1) = 228.7832;
params(22,1) = 228.7832;
params(23,1) = 45.7566; 
params(24,1) = 45.7566; 
params(25,1) = 2.4082;
params(26,1) = 2.4082;
params(27,1) = 1.2041e-08; 
params(28,1) = 1.0837e-06; 
params(28,2) = 2.4;
params(28,3) = 0.33;
params(28,4) = 45; 
params(29,1) = 0.0061; 
params(30,1) = 36.1237; 
params(30,2) = 1; 
params(31,1) = 0.0203; 
params(31,2) = 3.5; 
params(32,1) = 0.7225; 
params(32,2) = 3.5; 
params(33,1) = 0.1896; 
params(33,2) = 0.285714; 
params(34,1) = 0.0506; 
params(34,2) = 0.285714; 
params(35,1) = 0; 
params(35,2) = 3.5; 
params(36,1) = 0.4985; 
params(36,2) = 0.285714; 
params(37,1) = 0.0089;
params(38,1) = 0.0089; 
params(39,1) = 67.3467; 
params(40,1) = 67.3467; 
params(41,1) = 4.3059e-05; 
params(42,1) = 4.3059e-05; 
params(43,1) = 55.5617; 
params(44,1) = 55.5617;
params(45,1) = 45.7566; 
params(46,1) = 45.7566; 
params(47,1) = 1.6015; 
params(48,1) = 1.6015;
params(49,1) = 3.6124e-07; 
params(49,2) = 2.938; 
params(49,3) = 0.1775; 
params(49,4) = 45; 
params(50,1) = 0; 
params(50,2) = 3.5; 
params(51,1) = 0; 
params(51,2) = 3.5;
params(52,1) = 0.9970; 
params(52,2) = 0.285714; 
params(53,1) = 0.4985; 
params(53,2) = 0.285714; 
params(54,1) = 150.6280;
params(55,1) = 67.3467;
params(56,1) = 150.6280; 
params(57,1) = 67.3467;
params(58,1) = 0.0771;
params(59,1) = 0.0963; 
params(60,1) = 0.0770; 
params(61,1) = 0.0963;
params(62,1) = 228.7832;
params(63,1) = 55.5617; 
params(64,1) = 45.7566; 
params(65,1) = 45.7566; 
params(66,1) = 2.4082; 
params(67,1) = 1.6015; 

% TNF
params(68,1) = 0.0035; 
params(69,1) = 9.9027e-06;
params(70,1) = 0.0287; 
params(71,1) = 1.3245e3; 
params(71,2) = 0.001;
params(72,1) = 0.0253; 
params(72,2) = 0.001; 
params(73,1) = 0.1505;
params(74,1) = 41.0365;
params(75,1) = 0.0459; 
params(76,1) = 0.1505; 
params(77,1) = 1.6127;
params(78,1) = 385.6804; 
params(79,1) = 0.1505; 
params(80,1) = 2.2746e3; 
params(81,1) = 0.6247; 
params(82,1) = 22.6255; 
params(82,2) = 2;
params(82,3) = 0.001116; 

% IL1b
params(83,1) = 6.0206e-06; 
params(84,1) = 0.0181;
params(85,1) = 1.3245e3; 
params(85,2) = 0.001; 
params(86,1) = 0.0253; 
params(86,2) = 0.001; 
params(87,1) = 0.1204; 
params(88,1) = 48.1649;
params(89,1) = 0.0602; 
params(90,1) = 0.1204;
params(91,1) = 1.6127;
params(92,1) = 361.2367; 
params(93,1) = 0.1204;
params(94,1) = 2.047e3;



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


