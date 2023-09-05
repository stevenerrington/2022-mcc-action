%% Midcingulate analyses
% Extract sessions that only meet the following criteria:
dajo_datamap_curated = data_sessionCurate...
    (dajo_datamap,...
    'area', {'ACC'},... % - Sampled from midcingulate
    'monkey', {'dar','jou'},... % - In monkey Da and Jo
    'signal', {'SPK'},... % - Spike data only
    'spacing', [50, 100, 150]); % - From 50um, 100um, and 150um probes

% Find unique names for behavioral and neural files
dataFiles_beh = unique(dajo_datamap_curated.sessionBeh);
dataFiles_neural = unique(dajo_datamap_curated.dataFilename);

%% Behavioral analysis
% Extract behavior across all relevant sessions (as identified in the curated dajo_datamap).
% - Behavior includes: stop behavior, RTs, value behavior.
behavior = mcc_stopping_extractBeh(dirs,dataFiles_beh);

% Produce a summary figure
fig1_behavior_summary(dajo_datamap_curated, behavior, dataFiles_beh)

%% Neural analysis
% Define dorsal and ventral ACC boundaries
mcc_map_info = mcc_dv_mapCh(dajo_datamap_curated, behavior, dirs); % needs comments

% Extract spike density functions
% // mcc_getSpikeData

% In progress:
% // dev_initialPlot_stopping
% // dev_initialPlot_stoppingCluster
% // dev_cancel_glm

%% % % PLAN NOTES:

% (1)
% Single-neuron GLM. Adapt the analysis from Balewski et al., (2023)
% FR = β0 + β1 × trial type (no-stop or canceled) + β2 × SSD + β3 × trial number
% (https://github.com/t-elston/OFCvalue-to-ACCresponse/tree/main/analysis)
% aaa_Assess_OFC_DIRdecoding_UnitEncoding
% 100 ms bins, moved 10 ms; significance = p < 0.01 for at least 5
% consecutive bins (50 ms).

stopping_glm_analysis
ssd_glm_analysis

% (2)
% Cluster to discover common patterns of activity at the single-neuron
% level. Cluster solely on latency matched cancel/no-stop trials, not split
% by SSD

stopping_cluster_analysis

% (3)
% Population dynamics. Performed on significant neurons and entire
% sample.



% (4)
% Comparison between dMCC and vMCC

% (5)
% SSD "tuning" curves.

% (6)
% Dimensionality comparison between MCC and DMFC


