%% period_sensitivity.m
clear;
close all;
clc %clear session

%% Get current directory
fpath = mfilename('fullpath');
search = strfind(fpath, '/');
if isempty(search)
    search = strfind(fpath, '\');
end
fpath = fpath(1:search(end));
cd(fpath)
% Specify model directory
model_path = strcat(fpath(1:search(end-1)), "Temp_Model_Script/");
addpath(model_path)

% Baseline grouped parameters (1 x 62)
param_base = [ ...
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
    0.0058     % p68
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
    % IL1
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


ligand     = 'TNF';
% ligand     = 'IL1';

idx_groups = 1:numel(param_base);
delta      = 0.03;

% cd(fpath)
cd(model_path)

% Compute sensitivities
[c_vec, p0] = compute_period_sensitivities(param_base, ligand, idx_groups, delta);

cd(fpath)

% Save
fprintf('Baseline period p0 = %.2f hr\n', p0);
for k = 1:numel(idx_groups)
    fprintf('group %d: c_j = %.3f\n', idx_groups(k), c_vec(k));
end

c_TNF = c_vec;  
save('c_TNF.mat', 'c_TNF');
disp('Saved: c_TNF.mat');

% c_IL1 = c_vec;  
% save('c_IL1.mat', 'c_IL1');
% disp('Saved: c_IL1.mat');



%% =========================== Local functions ===========================

function [c_vec, p0] = compute_period_sensitivities(param_base, ligand, idx_groups, delta)
% compute_period_sensitivities
% param_base : 1 x 62 baseline grouped parameters
% ligand     : 'TNF' or 'IL1'
% idx_groups : indices of groups to compute sensitivities for
% delta      : relative perturbation size
    
    nG   = numel(idx_groups);
    c_vec = nan(nG,1);
    
    % baseline period
    p0 = get_period_from_params(param_base, ligand);

    for k = 1:nG
        j = idx_groups(k);
        kj0 = param_base(j);
        
        kj_up   = kj0 * (1 + delta);
        kj_down = kj0 * (1 - delta);

        par_up   = param_base;
        par_down = param_base;
        par_up(j)   = kj_up;
        par_down(j) = kj_down;

        % Periods
        p_up   = get_period_from_params(par_up,   ligand);
        p_down = get_period_from_params(par_down, ligand);

        % Log-based sensitivity: c_j ~= d log(P) / d log(k_j)
        num = log(p_up)   - log(p_down);
        den = log(kj_up)  - log(kj_down);
        c_vec(k) = num / den;
    end
end




function period_min = get_period_from_params(params, ligand)
% get_period_from_params
% params : 1 x 62 grouped parameters
% ligand : 'TNF' or 'IL1'

    % 1. Build p_mod
        p_mod = [
    1   1 params(1)   
    2   1 params(2)   
    3   1 params(3)   
    4   1 params(4)   
    5   1 params(5)   
    6   1 params(6)   
    7   1 params(7)   
    8   1 params(8)   
    30  1 params(8)   
    9   1 params(9)   
    10  1 params(10)  
    32  1 params(10)  
    11  1 params(11) 
    33  1 params(11)
    12  1 params(12)  
    34  1 params(12)
    13  1 params(13) 
    35  1 params(13)
    50  1 params(13)
    51  1 params(13)
    14  1 params(14)  
    52  1 params(14)
    15  1 params(15) 
    16  1 params(15)
    17  1 params(16)  
    18  1 params(16)
    19  1 params(17)  
    20  1 params(17)
    21  1 params(18) 
    22  1 params(18)
    62  1 params(18)
    23  1 params(19)  
    24  1 params(19)
    45  1 params(19)
    46  1 params(19)
    64  1 params(19)
    65  1 params(19)
    25  1 params(20)  
    26  1 params(20)
    66  1 params(20)
    27  1 params(21) 
    28  1 params(22)  
    29  1 params(23)  
    31  1 params(24)  
    36  1 params(25) 
    53  1 params(25)
    37  1 params(26)  
    38  1 params(26)
    39  1 params(27)  
    40  1 params(27)
    55  1 params(27)
    57  1 params(27)
    41  1 params(28)  
    42  1 params(28)
    43  1 params(29)  
    44  1 params(29)
    63  1 params(29)
    47  1 params(30) 
    48  1 params(30)
    67  1 params(30)
    49  1 params(31)  
    54  1 params(32)  
    56  1 params(32)
    58  1 params(33)  
    60  1 params(33)
    59  1 params(34)  
    61  1 params(34)

    % TNF
    68  1 params(35)
    69  1 params(36)
    70  1 params(37)
    71  1 params(38)
    72  1 params(39)
    73  1 params(40)
    74  1 params(41)
    75  1 params(42)
    76  1 params(43)
    77  1 params(44)
    78  1 params(45)
    79  1 params(46)
    80  1 params(47)
    81  1 params(48)
    82  1 params(49)

    % IL1
    83  1 params(50)
    84  1 params(51)
    85  1 params(52)
    86  1 params(53)
    87  1 params(54)
    88  1 params(55)
    89  1 params(56)
    90  1 params(57)
    91  1 params(58)
    92  1 params(59)
    93  1 params(60)
    94  1 params(61)
    ];


    % 2. Simulation settings
    names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA', ...
             'IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel', ...
             'IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln'};

    options = struct;
    options.DEBUG = 0;
    options.SIM_TIME = 12*60;
    time = linspace(0, options.SIM_TIME/60, options.SIM_TIME+1); 

    if strcmpi(ligand,'TNF')
        dose_scale = 1/5200;
    elseif strcmpi(ligand,'IL1')
        dose_scale = 1/17300;
    else
        error('ligand must be ''TNF'' or ''IL1''');
    end
    dose = 10;  
    
    [~, x] = temp_nfkbSimulate({ligand, dose*dose_scale}, names, p_mod, {}, options);
    
    rela_nuc = sum(x(:,6:8),2);
    
    [~, locs] = findpeaks(rela_nuc, time, 'MinPeakProminence', 0.01);
    if numel(locs) < 2
        period_min = NaN;
    else
        period_min = locs(2) - locs(1);  
    end
end



