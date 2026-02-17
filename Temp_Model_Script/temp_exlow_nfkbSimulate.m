function [t,x, simdata, x_end] = temp_exlow_nfkbSimulate(stimulus, names, params_mod, init_mod, options)
% Runs NF-kB model simulation with optional parameter and initial-state modifications.
% The function initializes the model, computes a basal steady state, then
% simulates the dynamic response to the specified stimulus.
% 24C

%% Default argument handling
if nargin <4
    init_mod = {};
    if nargin < 3
        params_mod = [];
        if nargin <2
            [~, tmp] = temp_exlow_nfkbInitialize;
            names = tmp.NAMES;
            if nargin<1
                stimulus = {};
            end
        end
    end
end

% OPTIONS
v.DT = 1;
v.SIM_TIME     = 18*60; 
v.T_EQUILIBRATE   = -4000;
v.DEBUG = 0;
v.PULSE_TIME = inf;
v.RELTOL = 1e-6;
v.ABSTOL = 1e-8;



if nargin>4
    fields = fieldnames(options);
    for i = 1:length(fields)
        v.(fields{i}) = options.(fields{i});
    end
end

if length(v.SIM_TIME) == 1
    v.SIM_TIME = [0 v.SIM_TIME];
end

if (v.SIM_TIME(1) > 0) && ~isfield(v,'STEADY_STATE')
    error('Initial/steady state vals (e.g. options.STEADY_STATE) must be provided if not starting simulation from t=0')
end


ode_opt = odeset('RelTol',v.RELTOL,'AbsTol',v.ABSTOL);


% Parameters and Initial values
if v.DEBUG;disp('INITIALIZATION: setting parameters and initial vals'); end
[v.PARAMS, v.SPECIES] = temp_exlow_nfkbInitialize();

for i = 1:size(params_mod,1)
    v.PARAMS(params_mod(i,1),params_mod(i,2)) = params_mod(i,3);
    if v.DEBUG
        disp(['- param(', num2str(params_mod(i,1)),',',num2str(params_mod(i,2)),...
            ') set to ',num2str(params_mod(i,3)) ])
    end
end

% Concentrations for fixed species 
v.SPECIES.INIT  = zeros(length(v.SPECIES.NAMES),1); % all species start at zero by default
v.SPECIES.INIT(strcmp('IKK_off',v.SPECIES.NAMES))   = 0.08;    
v.SPECIES.INIT(strcmp('RelA',v.SPECIES.NAMES))   = 0.08;    
v.SPECIES.INIT(strcmp('cRel',v.SPECIES.NAMES))   = 0.015;     
v.SPECIES.INIT(strcmp('TTR',v.SPECIES.NAMES))   = 0.035;    
v.SPECIES.INIT(strcmp('TAK1_off',v.SPECIES.NAMES))   = 0.1; 
v.SPECIES.INIT(strcmp('MIT',v.SPECIES.NAMES))   = 0.035;    


% Apply initial condition overrides
for i = 1:size(init_mod,1)
    idx = find(strcmp(init_mod{i,1},v.SPECIES.NAMES),1,'first');
    if ~isempty(idx)
    v.SPECIES.INIT(idx)   = init_mod{i,2};
    if v.DEBUG
        disp(['- initial ',init_mod{i,1},' altered to ',num2str(init_mod{i,2})])
    end
    elseif v.DEBUG
        disp(['- species ''',init_mod{i,1},''' not found']);
    end
end


%% Simulations

% PHASE1: find basal steady state
max_attempts = 100;
max_pct = 0.0005;
min_small = 8e-7; 

if ~isfield(v,'STEADY_STATE')
    v.PHASE     = 1;
    [~, r1] = ode15s('temp_nfkbOde', [v.T_EQUILIBRATE 0], v.SPECIES.INIT,ode_opt,v);
    starting_vals = r1(end,:);
    if v.DEBUG; disp('PHASE1: running model to basal steady state'); end
    converge = 0;
    iter = 0;
    temp_nfkbOde([],[],[],v); 
    while ~converge
        [~, r1] = ode15s('temp_nfkbOde', [v.T_EQUILIBRATE 0], starting_vals,ode_opt,v);
        pct_diff = abs((starting_vals - r1(end,:))./r1(end,:));
        starting_vals = r1(end,:);
        pct_diff(starting_vals < min_small) = 0; 
        if (max(pct_diff) < max_pct) || (iter>max_attempts)
            converge = 1;
            if (iter>max_attempts)
                idx = find(pct_diff==max(pct_diff),1,'first');
                disp(['Note: max steady state iterations reached (',num2str(max_attempts),'). Aborting with current state.'])
                disp(['[Non-convergent species: ', v.SPECIES.NAMES{idx},'(#',num2str(idx),')'...
                    '. Current state: ',num2str(starting_vals(idx)),' uM (',num2str(100*pct_diff(idx)),'% difference from prev.)]'])
            end
        end
        iter = iter+1;
    end
    v.STEADY_STATE = starting_vals;
    if v.DEBUG; disp(['- iterations required: ',num2str(iter),' (',num2str(-iter*v.T_EQUILIBRATE),' min)']); end
end

% PHASE2: add stimulus, and simulate from basal levels.
if v.DEBUG; disp('PHASE2: simulating dynamic response to stimulus'); end
starting_vals = v.STEADY_STATE;

for i = 1:size(stimulus,1)
    idx = find(strcmp(stimulus{i,1},v.SPECIES.NAMES),1,'first');
    if ~isempty(idx)
        starting_vals(idx)   = stimulus{i,2};
    if v.DEBUG
        disp(['- stimulus ''',stimulus{i,1},''' set to ',num2str(stimulus{i,2})])
    end
    elseif v.DEBUG
        disp(['- stimulus ',stimulus{i,1},' not found']);
    end
end

v.PHASE = 2;
[t_sim, x_sim] = ode15s('temp_nfkbOde', v.SIM_TIME,starting_vals,ode_opt,v);

if max(abs(imag(x_sim))) > eps
    warn(['Some numerical solutions have imaginary components. (Max imag. value = ',num2str(max(abs(imag(x_sim)))),'i)'])
end
x_sim = real(x_sim);
x_end = x_sim(end,:);


%% Interpolate outputs
x = [];
t = (v.SIM_TIME(1):v.DT:v.SIM_TIME(2))';

simdata = v;
simdata.STIMULUS = stimulus;

simdata.output = {};
for i = 1:length(names)
    idx = find(strcmp(names{i},v.SPECIES.NAMES),1,'first');
    if ~isempty(idx)
        x = cat(2,x,interp1(t_sim,x_sim(:,idx),t));
        simdata.output = cat(1,simdata.output,names{i});
    elseif v.DEBUG
        disp(['Output species ',names{i}, ' not found - skipping'])
    end
end

if v.DEBUG; disp('- - - - - - - - - - - - - - - - - - - -'); end


