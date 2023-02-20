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

sec1_figure_behaviorSummary



%% 2: 

