%% Get current directory
clear all; close all; clc %clear session

fpath = fileparts(mfilename('fullpath'));
cd(fpath);


% Specify Stimulations
run_TNF = true;
run_IL1 = true;


% Common settings
length_t = 12; % hr
names = {'IkBaRelA','IkBeRelA','IKKIkBaRelA','IKKIkBeRelA','RelA','IkBaRelAn','IkBeRelAn','RelAn','IkBacRel','IkBecRel','IKKIkBacRel','IKKIkBecRel','cRel','IkBacReln','IkBecReln','cReln'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = length_t*60;
time = linspace(0,options.SIM_TIME/60,options.SIM_TIME+1);
p_mod = [];


%% ---------------------------TNF SIMULATION ---------------------------
if run_TNF
    % TNF
    dose_scale_tnf = 1/5200; % Dose scaling from Shannon et al (2007)
    doses_tnf = [10];
    dose_tnf = doses_tnf(1) * dose_scale_tnf;

    % run 4 temps
    output_tnf_exlow = runOneTempSim(@temp_exlow_nfkbSimulate, 'TNF', dose_tnf, names, p_mod, options);
    output_tnf_low   = runOneTempSim(@temp_low_nfkbSimulate,   'TNF', dose_tnf, names, p_mod, options);
    output_tnf       = runOneTempSim(@temp_nfkbSimulate,       'TNF', dose_tnf, names, p_mod, options);
    output_tnf_high  = runOneTempSim(@temp_high_nfkbSimulate,  'TNF', dose_tnf, names, p_mod, options);

    % RelA Calculations (24C/34C/37C/40C)
    rela_nuc_tnf_exlow = sum(output_tnf_exlow(:,6:8),2);
    rela_nuc_tnf_low   = sum(output_tnf_low(:,6:8),2);
    rela_nuc_tnf       = sum(output_tnf(:,6:8),2);
    rela_nuc_tnf_high  = sum(output_tnf_high(:,6:8),2);

    % cRel Calculations (24C/34C/37C/40C)
    crel_nuc_tnf_exlow = sum(output_tnf_exlow(:,14:16),2);
    crel_nuc_tnf_low   = sum(output_tnf_low(:,14:16),2);
    crel_nuc_tnf       = sum(output_tnf(:,14:16),2);
    crel_nuc_tnf_high  = sum(output_tnf_high(:,14:16),2);

    % figure 1
    figure(1); clf
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % RelA plot
    nexttile; hold on
    plot(time, rela_nuc_tnf_exlow(:,1),'Color',[0 1 1],'LineWidth',2.5)
    plot(time, rela_nuc_tnf_low(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2.5)
    plot(time, rela_nuc_tnf(:,1),'Color',[0.9290 0.6940 0.1250], 'LineWidth',2.5)
    plot(time, rela_nuc_tnf_high(:,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5)
    axis([0 options.SIM_TIME/60 0 0.25]) % RelA
    xlabel('Time(hr)')
    ylabel('Nuclear Intensity')
    lgd = legend('24°C','34°C','37°C','40°C');
    fontsize(lgd,15,'points')
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.YLabel.FontSize = 15;
    ax.XAxis.FontSize = 15;

    % cRel plot
    nexttile; hold on
    plot(time, crel_nuc_tnf_exlow(:,1),'Color',[0 1 1],'LineWidth',2.5)
    plot(time, crel_nuc_tnf_low(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2.5)
    plot(time, crel_nuc_tnf(:,1),'Color',[0.9290 0.6940 0.1250], 'LineWidth',2.5)
    plot(time, crel_nuc_tnf_high(:,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5)
    axis([0 options.SIM_TIME/60 0 0.035]) % cRel
    xlabel('Time(hr)')
    ylabel('Nuclear Intensity')
    lgd = legend('24°C','34°C','37°C','40°C');
    fontsize(lgd,15,'points')
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.YLabel.FontSize = 15;
    ax.XAxis.FontSize = 15;
    set(gcf,'color','w');
end

%% ---------------------------IL1 SIMULATION ---------------------------
if run_IL1
    % IL1b
    doses_IL1 = [10]; % 10ng/mL
    dose_scale_IL1 = 1/17300; % IL1b molecular weight estimated between 17.3KDa and 31KDa
    dose_il1 = doses_IL1(1) * dose_scale_IL1;

    % run 4 temps
    output_IL1_exlow = runOneTempSim(@temp_exlow_nfkbSimulate, 'IL1', dose_il1, names, p_mod, options);
    output_IL1_low   = runOneTempSim(@temp_low_nfkbSimulate,   'IL1', dose_il1, names, p_mod, options);
    output_IL1       = runOneTempSim(@temp_nfkbSimulate,       'IL1', dose_il1, names, p_mod, options);
    output_IL1_high  = runOneTempSim(@temp_high_nfkbSimulate,  'IL1', dose_il1, names, p_mod, options);

    % RelA Calculations (24C/34C/37C/40C)
    rela_nuc_IL1_exlow = sum(output_IL1_exlow(:,6:8),2);
    rela_nuc_IL1_low   = sum(output_IL1_low(:,6:8),2);
    rela_nuc_IL1       = sum(output_IL1(:,6:8),2);
    rela_nuc_IL1_high  = sum(output_IL1_high(:,6:8),2);

    % cRel Calculations (24C/34C/37C/40C)
    crel_nuc_IL1_exlow = sum(output_IL1_exlow(:,14:16),2);
    crel_nuc_IL1_low   = sum(output_IL1_low(:,14:16),2);
    crel_nuc_IL1       = sum(output_IL1(:,14:16),2);
    crel_nuc_IL1_high  = sum(output_IL1_high(:,14:16),2);

    % figure 2
    figure(2); clf
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % RelA plot
    nexttile; hold on
    plot(time, rela_nuc_IL1_exlow(:,1),'Color',[0 1 1],'LineWidth',2.5)
    plot(time, rela_nuc_IL1_low(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2.5)
    plot(time, rela_nuc_IL1(:,1),'Color',[0.9290 0.6940 0.1250], 'LineWidth',2.5)
    plot(time, rela_nuc_IL1_high(:,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5)
    axis([0 options.SIM_TIME/60 0 0.25]) % RelA
    xlabel('Time(hr)')
    ylabel('Nuclear Intensity')
    lgd = legend('24°C','34°C','37°C','40°C');
    fontsize(lgd,15,'points')
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.YLabel.FontSize = 15;
    ax.XAxis.FontSize = 15;

    % cRel plot
    nexttile; hold on
    plot(time, crel_nuc_IL1_exlow(:,1),'Color',[0 1 1],'LineWidth',2.5)
    plot(time, crel_nuc_IL1_low(:,1),'Color',[0 0.4470 0.7410],'LineWidth',2.5)
    plot(time, crel_nuc_IL1(:,1),'Color',[0.9290 0.6940 0.1250], 'LineWidth',2.5)
    plot(time, crel_nuc_IL1_high(:,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5)
    axis([0 options.SIM_TIME/60 0 0.035]) % cRel
    xlabel('Time(hr)')
    ylabel('Nuclear Intensity')
    lgd = legend('24°C','34°C','37°C','40°C');
    fontsize(lgd,15,'points')
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.YLabel.FontSize = 15;
    ax.XAxis.FontSize = 15;
    set(gcf,'color','w');
end

%% =========================== Local function ===========================
function output = runOneTempSim(simFn, stimName, dose, names, p_mod, options)
% Minimal wrapper to avoid repeated blocks; keeps your original terms.
    [~, x] = simFn({stimName, dose}, names, p_mod, {}, options);
    output = x;
end