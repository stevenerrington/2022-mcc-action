
%% Setup workspace
clear all; clc; % Clear workspace
addpath('C:\Users\Steven\Desktop\Projects\2022-dajo-toolbox'); % DaJo-toolbox
addpath('C:\Users\Steven\Desktop\Projects\2022-acc-stopping\_toolbox\clustering_schall') % Clustering toolbox
dirs = data_setDir();
params.filter.band = [0.7 170]; params.filter.label = 'broadband';

%% Curate data
dajo_datamap = load_datamap(dirs);

dajo_datamap_curated = data_sessionCurate...
    (dajo_datamap,...
    'area', {'ACC'}, 'monkey', {'dar','jou'}, 'signal', {'SPK'}, 'spacing', [50, 100, 150]);

dataFiles_beh = unique(dajo_datamap_curated.sessionBeh);
dataFiles_neural = unique(dajo_datamap_curated.dataFilename);

%% Behavioral analysis
% Extract behavior
behavior = acc_stopping_extractBeh(dirs,dataFiles_beh);

% Look at RT adaptation following successful stopping
acc_stopping_RTadaptation


% Neural - spikes  -----------------------------------
signal_average_spk = acc_stopping_extractSDF(dirs,dataFiles_beh,dataFiles_neural,behavior);



%% Post-process neural signals
signal_collapse = neural_collapseSignalSession(signal_average_spk,...
    'events',{'target','stopSignal_artifical','ssrt'},...
    'conditions',{'C','GO'},...
    'conditions_map',[1 2]);

%% Clustering
acc_stopping_ssrtClustering





%% Archive
% %%
% % Neural - lfp     -----------------------------------
% acc_stopping_extractLFP(dirs,dataFiles_beh,dataFiles_neural,behavior,'filter',params.filter);
%         % 2022-05-20: This code takes a while and needs to be restructured.
