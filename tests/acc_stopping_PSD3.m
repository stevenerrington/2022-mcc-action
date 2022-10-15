%% Initialize
acc_depths = [-3000:50:3000];
n_freq = 129;
psd_array = nan(length(acc_depths),n_freq,n_sessions);

%% Analyize:
% Get session information
for session_i = 1:n_sessions
    acs_ch_mapping = util_getACCchannels(logInfo,session_i);
    neuralFilename = logInfo.neuralFilename{session_i};
    
    fprintf(['Extracting data for ' neuralFilename ': session %i of %i     \n'],...
        session_i, n_sessions);
    
    % Load LFP data and align on events
    data_in = load_lfpFile(dirs,neuralFilename);
    
    % Find faulty channels and extrapolate them to reduce contamation
    data_in = util_clearFaultyCh(data_in, neuralFilename, ephysLog);
    
    % Align LFP on target (3rd column in the trialEventTimes table).
    neurophys_lfp = lfp_alignTrials(behavior(beh_index).trialEventTimes(:,3), data_in, alignment_win);
    
    % Calculate PSD
    clear psd_analysis
    analysis_win = [-500:0];
    psd_analysis = lfp_getSessionPSD(neurophys_lfp, analysis_win+alignment_zero);
    
    
    for ch_i = 1:n_chan
        array_depth_index = find(acc_depths == acs_ch_mapping.um_depth_acs(ch_i));
        psd_array(array_depth_index,:,session_i) = psd_analysis.psd_norm(ch_i,:);
        
    end
end
