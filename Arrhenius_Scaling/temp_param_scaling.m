%% temp_param_scaling.m
% Arrhenius scaling of parameters across temperature and theoretical dp/dT.

clear; clc;

%% Inputs (edit these)
load('c_TNF.mat','c_TNF');   
load('c_IL1.mat','c_IL1');   


p_TNF37 = 110;   % period at 37C (min) / TNF - replace with your experimental measured value
p_IL137 = 131;   % period at 37C (min) / IL1 - replace with your experimental measured value

param37 = [ ...
    1          % p1
    2e-5       % p2
    2          % p3
    18         % p4
    5e-7       % p5
    5e-5       % p6
    0.0577623  % p7
    30         % p8
    0.0225     % p9
    0.6        % p10
    0.1575     % p11
    0.042      % p12
    0          % p13
    0.828      % p14
    0.0770164  % p15
    200        % p17
    8e-3       % p19
    190        % p21
    38         % p23
    2          % p25
    1e-8       % p27
    9e-07      % p28
    0.005064089% p29
    0.016875   % p31
    0.414      % p36
    0.00741282 % p37
    55.93008739% p39
    0.00003576 % p41
    46.1429    % p43
    1.33       % p47
    3e-7       % p49
    125.093633 % p54
    0.064      % p58
    0.08       % p59
    % TNF
    0.0029
%     0.0058     % p68
    8.224e-06  % p69
    0.02384    % p70
    1100       % p71
    0.021      % p72
    0.125      % p73
    34.08      % p74
    0.03812    % p75
    0.125      % p76
    1.3393     % p77
    320.3      % p78
    0.125      % p79
    1889       % p80
    0.5188     % p81
    18.79      % p82
    % IL-1
    5e-06      % p83
    0.015      % p84
    1100       % p85
    0.021      % p86
    0.1        % p87
    40         % p88
    0.05       % p89
    0.1        % p90
    1.3393     % p91
    300        % p92
    0.1        % p93
    1700       % p94
]';

m = numel(param37);

%% Example dp/dT from experiment
dpdT_exp_TNF_40_37  = -1;   % TNF 40–37°C
dpdT_exp_TNF_37_34  = -4;   % TNF 37–34°C

dpdT_exp_IL1_40_37  = -5;   % IL1 40–37°C
dpdT_exp_IL1_37_34  = -11;  % IL1 37–34°C



%  1. Arrhenius scaling template
R   = 8.3145e-3;           % kJ mol^-1 K^-1 
T37 = 310.15;              % 37°C in Kelvin
TvecC = [24 34 37 40];     % °C
TvecK = TvecC + 273.15;

E_vec = 50 * ones(m,1); % [kJ/mol]
A_vec = param37' .* exp(E_vec ./ (R*T37));
nT = numel(TvecK);
k_T = zeros(m, nT);  

for iT = 1:nT
    T = TvecK(iT);
    k_T(:,iT) = A_vec .* exp(-E_vec ./ (R*T));
end
% k_T(:,1) = 24°C, k_T(:,2) = 34°C, k_T(:,3) = 37°C, k_T(:,4) = 40°C



% 2. Theoretical dp/dT

T_eval = T37;

%% Remove NaN values from c_j vectors
valid_TNF = ~isnan(c_TNF);
valid_IL1 = ~isnan(c_IL1);

c_TNF_clean = c_TNF(valid_TNF);
c_IL1_clean = c_IL1(valid_IL1);

E_TNF_clean = E_vec(valid_TNF);
E_IL1_clean = E_vec(valid_IL1);

dpdT_theory = p_TNF37 * (1/(R*T_eval^2)) * (sum(c_TNF_clean .* E_TNF_clean)+sum(c_IL1_clean .* E_IL1_clean))/2;

% TNF
dpdT_TNF_theory = p_TNF37 * (1/(R*T_eval^2)) * sum(c_TNF_clean .* E_TNF_clean);
% IL1
dpdT_IL1_theory = p_IL137 * (1/(R*T_eval^2)) * sum(c_IL1_clean .* E_IL1_clean);




%  3. Save

save('Arrhenius_params_T.mat', 'k_T', 'TvecC', 'TvecK', ...
     'param37', 'A_vec', 'E_vec', ...
     'dpdT_TNF_theory', 'dpdT_IL1_theory');