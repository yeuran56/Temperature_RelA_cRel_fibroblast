%% Section 1: Declaration/initialization
function delta = temp_nfkbOde(t,x,ode_options,v)

% Initialize memory used for delayed reactions in the ODE model
persistent DELAY;

if isempty(t)
    sz = 5000;  % Size of total delay memory
    DELAY.t = zeros(sz,1);
    DELAY.rela = zeros(sz,1); 
    DELAY.ikbat = zeros(sz,1); 
    DELAY.crel = zeros(sz,1); 
    DELAY.ikbet = zeros(sz,1); 
    DELAY.idx     = 1; 
    return;
end


% Slice parameters, get previous concentrations
p = v.PARAMS;
delta = zeros(size(x));



%% Section 2: Unpack species
% TNF and IL1
stim = x(1);
IkBa = x(2);
IkBan = x(3);
IkBe = x(4);
IkBen = x(5);
IkBaRelA = x(6);
IkBaRelAn = x(7);
IkBacRel = x(8);
IkBacReln = x(9);
IkBeRelA = x(10);
IkBeRelAn = x(11);
IkBecRel = x(12);
IkBecReln = x(13);
IkBat = x(14);
IkBet = x(15);
IKKIkBaRelA = x(16);
IKKIkBacRel = x(17);
IKKIkBa = x(18);
IKKIkBeRelA = x(19);
IKKIkBecRel = x(20);
IKKIkBe = x(21);
cRel = x(22);
cReln = x(23);
RelA = x(24);
RelAn = x(25);
IKK_off = x(26);
IKK = x(27);
IKK_i = x(28);

% TNF
TNF = x(29);
TNFR = x(30);
TNFR_TNF = x(31);
TTR = x(32);
C1_off = x(33);
C1 = x(34);
TAK1_off = x(35);
TAK1 = x(36);

% IL1
IL1 = x(37);
IL1R = x(38);
IL1R_IL1 = x(39);
MIT = x(40);
C2_off = x(41);
C2 = x(42);




%% Section 3: Calculate delays
% --- Discrete delay handling ------------------------------------------------
% This section implements fixed time delays for selected reactions by
% caching past values and interpolating x(t - tau).
% Note: delay targets/tau indices are hard-coded; update this block if the
% reaction ordering or parameter indexing changes.

RelAn_tau6   = RelAn;
IkBat_tau8   = IkBat;
cReln_tau28   = cReln;
IkBet_tau30   = IkBet;
RelAn_tau49   = RelAn;
if v.PHASE ~= 1
    idx = find(DELAY.t(1:DELAY.idx)>=t,1,'first');
    if isempty(idx)
        idx = DELAY.idx+1;
    end
    DELAY.t(idx) = t;
    DELAY.idx = idx;
    DELAY.rela(DELAY.idx) = RelAn;
    DELAY.ikbat(DELAY.idx) = IkBat;
    DELAY.crel(DELAY.idx) = cReln;
    DELAY.ikbet(DELAY.idx) = IkBet;
    if DELAY.idx > length(DELAY.t)
        error(['Delay index out of range: increase memory size (line 13 - ',...
            'sz = ',num2str(length(DELAY.t)),';). Reached t = ',num2str(DELAY.t(end))])
    end

    % 1. IkBa mRNA transcription/processing delay (RXN 6)
    tau = p(6,4);
    if tau > 0 
        if t > tau
            RelAn_tau6 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.rela(1:DELAY.idx-1),t-tau);
        else
            RelAn_tau6 = ...
                DELAY.rela(1); 
        end
    end
    
    % 2. IkBa protein translation/maturation delay (RXN 8)
    tau = p(8,2);
    if tau > 0 
        if t > tau
            IkBat_tau8 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.ikbat(1:DELAY.idx-1),t-tau);
        else
            IkBat_tau8 = ...
                DELAY.ikbat(1); 
        end
    end

    % 3. IkBe mRNA transcription/processing delay (RXN 28)
    tau = p(28,4);
    if tau > 0  
        if t > tau
            cReln_tau28 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.crel(1:DELAY.idx-1),t-tau);
        else
            cReln_tau28 = ...
                DELAY.crel(1); 
        end
    end

    % 4. IkBe protein translation/maturation delay (RXN 8)
    tau = p(30,2);
    if tau > 0  
        if t > tau
            IkBet_tau30 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.ikbet(1:DELAY.idx-1),t-tau);
        else
            IkBet_tau30 = ...
                DELAY.ikbet(1);
        end
    end
    
    % 5. IkBa mRNA transcription/processing delay (RXN 49)
    tau = p(49,4);
    if tau > 0  
        if t > tau
            RelAn_tau49 = ...
                interp1(DELAY.t(1:DELAY.idx-1),...
                        DELAY.rela(1:DELAY.idx-1),t-tau);
        else
            RelAn_tau49 = ...
                DELAY.rela(1);
        end
    end
  
end

% After the pulse duration, external stimuli are set to zero
if isfield(v,'PULSE_TIME') && (t > v.PULSE_TIME)
    TNF = 0;
    IL1 = 0;
    stim=0;
end



%% Section 4: Set reaction rates
rxn_1 = p(1,1) * stim * IKK_off;
rxn_2 = p(2,1) * IKK_off;
rxn_3 = p(3,1) * IKK;
rxn_4 = p(4,1) * IKK_i;
rxn_5 = p(5,1);
rxn_6 = p(6,1) * (RelAn_tau6.^p(6,2))/( (RelAn_tau6.^p(6,2)) + (p(6,3).^p(6,2)) );
rxn_7 = p(7,1) * IkBat;
rxn_8 = p(8,1) * IkBat_tau8;
rxn_9 = p(9,1) * IkBa;
rxn_10 = p(10,1) * RelA;
rxn_11 = p(11,1) * IkBan;
rxn_12 = p(12,1) * RelAn;
rxn_13 = p(13,1) * IkBaRelA;
rxn_14 = p(14,1) * IkBaRelAn;
rxn_15 = p(15,1) * IkBa;
rxn_16 = p(16,1) * IkBan;
rxn_17 = p(17,1) * IkBa * RelA;
rxn_18 = p(18,1) * IkBan * RelAn;
rxn_19 = p(19,1) * IkBaRelA;
rxn_20 = p(20,1) * IkBaRelAn;
rxn_21 = p(21,1) * IKK * IkBaRelA;
rxn_22 = p(22,1) * IKK * IkBa;
rxn_23 = p(23,1) * IKKIkBaRelA;
rxn_24 = p(24,1) * IKKIkBa;
rxn_25 = p(25,1) * IKKIkBaRelA;
rxn_26 = p(26,1) * IKKIkBa;
rxn_27 = p(27,1);
rxn_28 = p(28,1) * (cReln_tau28.^p(28,2))/( (cReln_tau28.^p(28,2)) + (p(28,3).^p(28,2)) );
rxn_29 = p(29,1) * IkBet;
rxn_30 = p(30,1) * IkBet_tau30;
rxn_31 = p(31,1) * IkBe;
rxn_32 = p(32,1) * cRel;
rxn_33 = p(33,1) * IkBen;
rxn_34 = p(34,1) * cReln;
rxn_35 = p(35,1) * IkBecRel;
rxn_36 = p(36,1) * IkBecReln;
rxn_37 = p(37,1) * IkBe;
rxn_38 = p(38,1) * IkBen;
rxn_39 = p(39,1) * IkBe * cRel;
rxn_40 = p(40,1) * IkBen * cReln;
rxn_41 = p(41,1) * IkBecRel;
rxn_42 = p(42,1) * IkBecReln;
rxn_43 = p(43,1) * IKK * IkBecRel;
rxn_44 = p(44,1) * IKK * IkBe;
rxn_45 = p(45,1) * IKKIkBecRel;
rxn_46 = p(46,1) * IKKIkBe;
rxn_47 = p(47,1) * IKKIkBecRel;
rxn_48 = p(48,1) * IKKIkBe;
rxn_49 = p(49,1) * (RelAn_tau49.^p(49,2))/( (RelAn_tau49.^p(49,2)) + (p(49,3).^p(49,2)) );
rxn_50 = p(50,1) * IkBacRel;
rxn_51 = p(51,1) * IkBeRelA;
rxn_52 = p(52,1) * IkBacReln;
rxn_53 = p(53,1) * IkBeRelAn;
rxn_54 = p(54,1) * IkBa * cRel;
rxn_55 = p(55,1) * IkBe * RelA;
rxn_56 = p(56,1) * IkBan * cReln;
rxn_57 = p(57,1) * IkBen * RelAn;
rxn_58 = p(58,1) * IkBacRel;
rxn_59 = p(59,1) * IkBeRelA;
rxn_60 = p(60,1) * IkBacReln;
rxn_61 = p(61,1) * IkBeRelAn;
rxn_62 = p(62,1) * IKK * IkBacRel;
rxn_63 = p(63,1) * IKK * IkBeRelA;
rxn_64 = p(64,1) * IKKIkBacRel;
rxn_65 = p(65,1) * IKKIkBeRelA;
rxn_66 = p(66,1) * IKKIkBacRel;
rxn_67 = p(67,1) * IKKIkBeRelA;

% TNF
rxn_68 = p(68,1) * TNF;
rxn_69 = p(69,1);
rxn_70 = p(70,1) * TNFR;
rxn_71 = p(71,1) * TNF * TNFR;
rxn_72 = p(72,1) * TNFR_TNF;
rxn_73 = p(73,1) * TNFR_TNF;
rxn_74 = p(74,1) * TNFR_TNF * TTR;
rxn_75 = p(75,1) * C1_off;
rxn_76 = p(76,1) * C1_off;
rxn_77 = p(77,1) * C1_off;
rxn_78 = p(78,1) * C1;
rxn_79 = p(79,1) * C1;
rxn_80 = p(80,1) * C1 * TAK1_off;
rxn_81 = p(81,1) * TAK1;
rxn_82 = p(82,1) * (TAK1.^p(82,2))/( (TAK1.^p(82,2)) + (p(82,3).^p(82,2)) ) * IKK_off;

% IL1
rxn_83 = p(83,1);
rxn_84 = p(84,1) * IL1R;
rxn_85 = p(85,1) * IL1 * IL1R;
rxn_86 = p(86,1) * IL1R_IL1;
rxn_87 = p(87,1) * IL1R_IL1;
rxn_88 = p(88,1) * IL1R_IL1 * MIT;
rxn_89 = p(89,1) * C2_off;
rxn_90 = p(90,1) * C2_off;
rxn_91 = p(91,1) * C2_off;
rxn_92 = p(92,1) * C2;
rxn_93 = p(93,1) * C2;
rxn_94 = p(94,1) * C2 * TAK1_off;



%% Section 5: Set species' deltas from reactions
delta(1) = 0;
delta(2) = - rxn_9 - rxn_15 - rxn_17 - rxn_22 - rxn_54 + rxn_8 + (rxn_11 * p(11,2)) + rxn_19 + rxn_24 + rxn_58;
delta(3) = - rxn_11 - rxn_16 - rxn_18 - rxn_56 + (rxn_9 * p(9,2)) + rxn_20 + rxn_60;
delta(4) = - rxn_31 - rxn_37 - rxn_39 - rxn_44 - rxn_55 + rxn_30 + (rxn_33 * p(33,2)) + rxn_41 + rxn_46 + rxn_59;
delta(5) = - rxn_33 - rxn_38 - rxn_40 - rxn_57 + (rxn_31 * p(31,2)) + rxn_42 + rxn_61;
delta(6) = - rxn_13 - rxn_19 - rxn_21 + (rxn_14 * p(14,2)) + rxn_17 + rxn_23;
delta(7) = - rxn_14 - rxn_20 + (rxn_13 * p(13,2)) + rxn_18;
delta(8) = - rxn_50 - rxn_58 - rxn_62 + (rxn_52 * p(52,2)) + rxn_54 + rxn_64;
delta(9) = - rxn_52 - rxn_60 + (rxn_50 * p(50,2)) + rxn_56;
delta(10) = - rxn_51 - rxn_59 - rxn_63 + (rxn_53 * p(53,2)) + rxn_55 + rxn_65;
delta(11) = - rxn_53 - rxn_61 + (rxn_51 * p(51,2)) + rxn_57;
delta(12) = - rxn_35 - rxn_41 - rxn_43 + (rxn_36 * p(36,2)) + rxn_39 + rxn_45;
delta(13) = - rxn_36 - rxn_42 + (rxn_35 * p(35,2)) + rxn_40;
delta(14) = - rxn_7 + rxn_5 + rxn_6;
delta(15) = - rxn_29 + rxn_27 + rxn_28 + rxn_49;
delta(16) = - rxn_23 - rxn_25 + rxn_21;
delta(17) = - rxn_64 - rxn_66 + rxn_62;
delta(18) = - rxn_24 - rxn_26 + rxn_22;
delta(19) = - rxn_65 - rxn_67 + rxn_63;
delta(20) = - rxn_45 - rxn_47 + rxn_43;
delta(21) = - rxn_46 - rxn_48 + rxn_44;
delta(22) = - rxn_32 - rxn_39 - rxn_54 + (rxn_34 * p(34,2)) + rxn_41 + rxn_47 + rxn_58 + rxn_66;
delta(23) = - rxn_34 - rxn_40 - rxn_56 + (rxn_32 * p(32,2)) + rxn_42 + rxn_60;
delta(24) = - rxn_10 - rxn_17 - rxn_55 + (rxn_12 * p(12,2)) + rxn_19 + rxn_25 + rxn_59 + rxn_67;
delta(25) = - rxn_12 - rxn_18 - rxn_57 + (rxn_10 * p(10,2)) + rxn_20 + rxn_61;
delta(26) = - rxn_1 - rxn_2 - rxn_82 + rxn_4;
delta(27) = - rxn_3 - rxn_21 - rxn_22 - rxn_43 - rxn_44 - rxn_62 - rxn_63 + rxn_1 + rxn_2 + rxn_23 + rxn_24 + rxn_25 + rxn_26 + rxn_45 + rxn_46 + rxn_47 + rxn_48 + rxn_64 + rxn_65 + rxn_66 + rxn_67 + rxn_82;
delta(28) = - rxn_4 + rxn_3;

% TNF
delta(29) = - rxn_68 - (rxn_71 * p(71,2)) + (rxn_72 * p(72,2));
delta(30) = - rxn_70 - rxn_71 + rxn_69 + rxn_72;
delta(31) = - rxn_72 - rxn_73 - rxn_74 + rxn_71 + rxn_75;
delta(32) = - rxn_74 + rxn_75;
delta(33) = - rxn_75 - rxn_76 - rxn_77 + rxn_74 + rxn_78;
delta(34) = - rxn_78 - rxn_79 + rxn_77;
delta(35) = - rxn_80 + rxn_81 - rxn_94;
delta(36) = - rxn_81 + rxn_80 + rxn_94;

% IL1b
delta(37) = - (rxn_85 * p(85,2)) + (rxn_86 * p(86,2));
delta(38) = - rxn_84 - rxn_85 + rxn_83 + rxn_86;
delta(39) = - rxn_86 - rxn_87 - rxn_88 + rxn_85 + rxn_89;
delta(40) = - rxn_88 + rxn_89;
delta(41) = - rxn_89 - rxn_90 - rxn_91 + rxn_88 + rxn_92;
delta(42) = - rxn_92 - rxn_93 + rxn_91;

end
