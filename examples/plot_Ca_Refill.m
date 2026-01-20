%% Compare reservoir refill traces: tubular (diamond) vs bubble (HCP) for R10
% This script:
%   1) Loads ONE mesh file (any replicate is fine) to provide mesh structure
%   2) Parses snapshot files using parseReservoirSims
%   3) Saves compact .mat outputs for reuse
%   4) Loads those .mat files and plots ER free Ca (normalized) vs time

clear; clc;

%% -------------------- USER SETTINGS --------------------
% Tubular (diamond) directory + file naming
dir_tub   = 'example/';
mesh_tub  = 'diamondR10_nf_nuccon91.0.0.mesh.txt';
glob_tub  = 'diamondR10_nf_nuccon91.%d.%d.snap.txt';
out_tub   = 'diamondR10_test.mat';

% Bubble (HCP) directory + file naming
dir_bub   = 'example/';
mesh_bub  = 'HCP_nf_bub_10sheet_buffscl_R10.0.0.mesh.txt';
glob_bub  = 'HCP_nf_bub_10sheet_buffscl_R10.%d.%d.snap.txt';
out_bub   = 'HCP_nf_bub_10sheet_buffscl_R10_test.mat';

% Snapshot parsing (here only 1 param, 1 replicate — adjust later if needed)
nparam = 1;
nrep   = 1;

% Scaling to convert the stored 1D quantity back into "relative to max"
% (you were using 1/(0.5*pi*r^2) for both)
scl_tub = 1/(0.5*pi*0.05^2);    % tubular radius = 0.05 um
scl_bub = 1/(0.5*pi*0.018^2);   % bubble tube radius = 0.018 um

%% -------------------- STEP 1: PARSE + SAVE (TUBULAR) --------------------
% Load one mesh (any run is fine; mesh geometry should match the setup)
MSHs = MeshObj(fullfile(dir_tub, mesh_tub));   % a single MeshObj is OK
fileglob = fullfile(dir_tub, glob_tub);

[rvalavg,snaptimes,rvalall,avgfieldsavg,avgfieldsall] = ...
    parseReservoirSims(fileglob, nparam, nrep, MSHs);

save(fullfile(dir_tub, out_tub), ...
    'rvalavg','snaptimes','rvalall','nparam','nrep','avgfieldsavg','avgfieldsall');

%% -------------------- STEP 2: PARSE + SAVE (BUBBLE) --------------------
MSHs = MeshObj(fullfile(dir_bub, mesh_bub));
fileglob = fullfile(dir_bub, glob_bub);

[rvalavg,snaptimes,rvalall,avgfieldsavg,avgfieldsall] = ...
    parseReservoirSims(fileglob, nparam, nrep, MSHs);

save(fullfile(dir_bub, out_bub), ...
    'rvalavg','snaptimes','rvalall','nparam','nrep','avgfieldsavg','avgfieldsall');

%% -------------------- STEP 3: LOAD SAVED RESULTS --------------------
D_tub = load(fullfile(dir_tub, out_tub));
t_tub = D_tub.snaptimes;
u_tub = D_tub.avgfieldsavg;   % assumes 1D time series already

D_bub = load(fullfile(dir_bub, out_bub));
t_bub = D_bub.snaptimes;
u_bub = D_bub.avgfieldsavg;

%% -------------------- STEP 4: PLOT --------------------
figure; hold on;

% Tubular: solid
plot(t_tub, u_tub * scl_tub, ...
    'LineWidth', 2, ...
    'LineStyle', '-', ...
    'DisplayName', 'tubular (diamond)');

% Bubble: dashed
plot(t_bub, u_bub * scl_bub, ...
    'LineWidth', 2, ...
    'LineStyle', '--', ...
    'DisplayName', 'bubble (HCP)');

plot_cleanup(gca,'FontSize',14);
xlabel('time (sec)');
ylabel('ER free Ca$^{2+}$, relative to max');
ylim([0 1]);
axis square;
legend('Location','best');
title('Reservoir refill comparison (R10)');

