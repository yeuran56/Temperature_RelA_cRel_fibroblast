function [params, species] = temp_nfkbInitialize()
% Initialize parameter set and species definitions for the NF-kB ODE model.
% 37C

% Core pathway parameters
params(1,1) = 1; % induced IKK activation
params(2,1) = 2e-05; % basal IKK activation
params(3,1) = 2; % IKK inhibition
params(4,1) = 18; % IKK renewal 
params(5,1) = 5e-7; % basal IkBa mRNA synthesis
params(6,1) = 5e-05; % induced IkBa mRNA synthesis rate
params(6,2) = 2.4; % Hill coefficient
params(6,3) = 0.33; % EC50
params(6,4) = 14; % transcription/processing delay
params(7,1) = 0.0577623; % IkBa mRNA degradation
params(8,1) = 30; % IkBa mRNA translation
params(8,2) = 1; % Translation/folding delay
params(9,1) = 0.0225; % nuclear import of IkBa
params(9,2) = 3.5; % Volume scale
params(10,1) = 0.6; % nuclear import of RelA
params(10,2) = 3.5; % Volume scale
params(11,1) = 0.1575; % nuclear export of IkBa
params(11,2) = 0.285714; % Volume scale
params(12,1) = 0.042; % nuclear export of RelA
params(12,2) = 0.285714; % Volume scale
params(13,1) = 0; % nuclear import of IkBa-RelA
params(13,2) = 3.5; % Volume scale
params(14,1) = 0.828; % nuclear export of IkBa-RelA
params(14,2) = 0.285714; % Volume scale
params(15,1) = 0.0770164; % degradation of IkBa
params(16,1) = 0.0770164; % degradation of IkBa (nuc)
params(17,1) = 200; % IkBa-RelA association
params(18,1) = 200; % IkBa-RelA association (nuc)
params(19,1) = 0.008; % IkBa-RelA dissociation
params(20,1) = 0.008; % IkBa-RelA dissociation (nuc)
params(21,1) = 190; % IKK-IkBa-RelA association
params(22,1) = 190; % IKK-IkBa association
params(23,1) = 38; % IKK-IkBa-RelA dissociation
params(24,1) = 38; % IKK-IkBa dissociation
params(25,1) = 2; % Phosphorylation/degradation of complexed IkBa
params(26,1) = 2; % Phosphorylation/degradation of IkBa
params(27,1) = 1e-08; % basal IkBe mRNA synthesis
params(28,1) = 9e-07; % induced IkBe mRNA synthesis from cRel
params(28,2) = 2.4; % Hill coefficient
params(28,3) = 0.33; % EC50
params(28,4) = 45; % mRNA transcription/processing/maturation delay
params(29,1) = 0.00506409; % IkBe mRNA degradation
params(30,1) = 30; % IkBe mRNA translation
params(30,2) = 1; % Translation/folding delay
params(31,1) = 0.016875; % nuclear import of IkBe
params(31,2) = 3.5; % Volume scale
params(32,1) = 0.6; % nuclear import of cRel:p50
params(32,2) = 3.5; % Volume scale
params(33,1) = 0.1575; % nuclear export of IkBe
params(33,2) = 0.285714; % Volume scale
params(34,1) = 0.042; % nuclear export of cRel:p50
params(34,2) = 0.285714; % Volume scale
params(35,1) = 0; % nuclear import of IkBe-NFkB
params(35,2) = 3.5; % Volume scale
params(36,1) = 0.414; % nuclear export of IkBe-NFkB
params(36,2) = 0.285714; % Volume scale
params(37,1) = 0.00741282; % degradation of IkBe
params(38,1) = 0.00741282; % degradation of IkBe (nuc)
params(39,1) = 55.9301; % IkBe-cRel association
params(40,1) = 55.9301; % IkBe-cRel association (nuc)
params(41,1) = 3.576e-05; % IkBe-cRel dissociation
params(42,1) = 3.576e-05; % IkBe-cRel dissociation (nuc)
params(43,1) = 46.1429; % IKK-IkBe-cRel association
params(44,1) = 46.1429; % IKK-IkBe association
params(45,1) = 38; % IKK-IkBe-cRel dissociation
params(46,1) = 38; % IKK-IkBe dissociation
params(47,1) = 1.33; % Phosphorylation/degradation of complexed IkBe
params(48,1) = 1.33; % Phosphorylation/degradation of IkBe
params(49,1) = 3e-07; % induced IkBe mRNA synthesis from RelA
params(49,2) = 2.938; % Hill coefficient
params(49,3) = 0.1775; % EC50 
params(49,4) = 45; % mRNA transcription/processing/maturation delay
params(50,1) = 0; % nuclear import of IkBa-NFkB
params(50,2) = 3.5; % Volume scale
params(51,1) = 0; % nuclear import of IkBe-NFkB
params(51,2) = 3.5; % Volume scale
params(52,1) = 0.828; % nuclear export of IkBa-NFkB
params(52,2) = 0.285714; % Volume scale
params(53,1) = 0.414; % nuclear export of IkBe-NFkB
params(53,2) = 0.285714; % Volume scale
params(54,1) = 125.094; % IkBa-cRel association
params(55,1) = 55.9301; % IkBe-RelA association
params(56,1) = 125.094; % IkBa-cRel association (nuc)
params(57,1) = 55.9301; % IkBe-RelA association (nuc)
params(58,1) = 0.064; % IkBa-cRel dissociation
params(59,1) = 0.08; % IkBe-RelA dissociation
params(60,1) = 0.064; % IkBa-cRel dissociation (nuc)
params(61,1) = 0.08; % IkBe-RelA dissociation (nuc)
params(62,1) = 190; % IKK-IkBa-cRel association
params(63,1) = 46.1429; % IKK-IkBe-RelA association
params(64,1) = 38; % IKK-IkBa-cRel dissociation
params(65,1) = 38; % IKK-IkBe-RelA dissociation
params(66,1) = 2; % Phosphorylation/degradation of complexed IkBa
params(67,1) = 1.33; % Phosphorylation/degradation of complexed IkBe

% TNF
params(68,1) = 0.0029; % TNF degradation
params(69,1) = 8.224e-06; % TNFR synthesis
params(70,1) = 0.02384; % TNFR degradation
params(71,1) = 1100; % Capture of TNF
params(71,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(72,1) = 0.021; % Release of TNF
params(72,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(73,1) = 0.125; % Internalization/degradation of complexed TNFR
params(74,1) = 34.08; % Association of TRAF2/RIP1 with receptor trimer
params(75,1) = 0.03812; % Dissociation of TRAF2/RIP1 with receptor trimer
params(76,1) = 0.125; % Internalization/degradation of complexed TNFR 
params(77,1) = 1.3393; % Activation of complexed TNFR
params(78,1) = 320.3; % Inactivation of complexed TNFR
params(79,1) = 0.125; % Internalization/degradation of complexed TNFR
params(80,1) = 1889; % Activation of TAK1 by C1
params(81,1) = 0.5188; % Inactivation of TAK1
params(82,1) = 18.79; % Activation of IKK by TAK1
params(82,2) = 2; % Hill coefficient for IKK activation 
params(82,3) = 0.001116; % EC50 for mRNA syn

% IL1b
params(83,1) = 5e-06; % IL1R synthesis
params(84,1) = 0.015; % IL1R degradation
params(85,1) = 1100; % Capture of IL1
params(85,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(86,1) = 0.021; % Release of IL1 
params(86,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(87,1) = 0.1; % Internalization/removal rate of ligand-bound IL1R complex
params(88,1) = 40; % Recruitment rate of adaptor module (MIT) to activated IL1R
params(89,1) = 0.05; % Dissociation rate of adaptor module (MIT) from receptor complex
params(90,1) = 0.1; % Turnover of adaptor-bound IL1R signaling complex
params(91,1) = 1.3393; % Activation of complexed IL1R
params(92,1) = 300; % Inactivation of complexed IL1R
params(93,1) = 0.1; % Internalization/degradation of complexed IL1R
params(94,1) = 1700; % Activation of TAK1 by C2


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