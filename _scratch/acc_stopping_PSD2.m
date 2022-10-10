
alignment_win = [-1000 2000]; alignment_zero = abs(alignment_win(1));
analysis_win = [-500:0];

session_i = 1;

% Find session details
neuralFilename = dajo_datamap_curated.dataFilename{session_i};
behFilename = data_findBehFile(neuralFilename);
beh_index = util_find_beh_index(behavior,behFilename);
fprintf(['Analysing session %i of %i: ' neuralFilename '    \n'], session_i, n_sessions)

% Load LFP data and align on events
data_in = load_lfpFile(dirs,neuralFilename);
neurophys_lfp = lfp_alignTrials(behavior(beh_index).trialEventTimes(:,3), data_in, alignment_win);

% Calculate PSD
psd_analysis{session_i,1} = lfp_getSessionPSD(neurophys_lfp, analysis_win+alignment_zero);

