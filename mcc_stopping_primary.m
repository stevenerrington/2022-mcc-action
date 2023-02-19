%{
Properties of midcingulate neurons in a saccade countermanding task
2022-10-11

Dependencies:
- gramm toolbox
%}

%% Setup workspace
clear all; clc; % Clear workspace
dirs = get_dirs_mcc('home'); 

getColors_mcc

% Setup directories and run global functions
ephysLog = importOnlineEphysLogMaster;

%% Curate data
% Load in the main pre-processed data map
dajo_datamap = load_datamap(fullfile(dirs.root,'data'));

% Curate the datamap to only have sessions of interest
%   Criteria: ACC sessions, in Da and Jo, that contain spikes, across all
%   all electrode spacings
dajo_datamap_curated = data_sessionCurate...
    (dajo_datamap,...
    'area', {'ACC'}, 'monkey', {'dar','jou'}, 'signal', {'SPK'}, 'spacing', [50, 100, 150]);

% Once curated, get the names of the relevant behavior and neural files
dataFiles_beh = unique(dajo_datamap_curated.sessionBeh);
dataFiles_neural = unique(dajo_datamap_curated.dataFilename);

%% 1: Behavioral analysis
% Extract behavior across all relevant sessions (as identified in the
% curated dajo_datamap).
%   Behavior includes: stop behavior, RTs, value behavior.
behavior = mcc_stopping_extractBeh(dirs,dataFiles_beh);

% Look at RT adaptation following successful stopping
%   We approach this in two ways:
%       1. The median RT on no-stop trials following canceled, no-stop, and
%       non-canceled trials
%       2. The change in RT from a no-stop trial prior to a canceled trial,
%       to a no-stop trial following a canceled trial.
acc_beh_RTadaptation


%% 2: Dorsal/ventral banks
acc_dv_mapCh
% acc_getSpikeData

%% 3: Functionality (action)
% 3.1 Error
% 3.1.1 Extract saccade-aligned spike density functions


% 3.1.2 Find neurons that 






% 3.2 Stopping
% 3.2.1 Extract stopping spike density functions

% 3.1.2 Get average SDF across different SSDs

% 3.1.3 Run clustering algorithm
%   This clusters on activity aligned on stop-signal (0 to 200 ms)






%% 4: Functionality (value)



