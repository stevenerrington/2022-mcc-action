function behavior = acc_stopping_extractBeh(dirs,dataFiles_beh)

% Looping through each of the individual data files
parfor dataFileIdx = 1:length(dataFiles_beh)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_beh{dataFileIdx},dataFileIdx,length(dataFiles_beh))
    
    % We then get the (behavior) filename of the record of interest
    behFilename = [dataFiles_beh{dataFileIdx} '-beh'];
    % and load it into the workspace
    import_data = struct(); import_data = load_behFile(dirs,behFilename);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Once we have the behavior in matlab, we can then look to extract
    % relevant task/session behavior.
    
    sessionName = {behFilename}; % Get the session name
    [ttx, ttx_history, trialEventTimes] =... % Index of trials and event timings
        beh_getTrials(import_data.events.stateFlags_,import_data.events.Infos_);
    [stopSignalBeh, ~] = beh_getStoppingInfo... % Stopping behavior
        (import_data.events.stateFlags_,import_data.events.Infos_,ttx);
    trialEventTimes.ssrt = trialEventTimes.stopSignal_artifical + stopSignalBeh.ssrt.integrationWeighted;
    [ttm] = beh_getMatchedTrials(stopSignalBeh,ttx, trialEventTimes); % Trial matching indices
        
    % After extracting the individual behavioral variable, we then collapse
    % it into one structure for the given session.
    behavior(dataFileIdx) = struct('sessionName',sessionName,'ttx',ttx,'trialEventTimes',trialEventTimes,...
        'stopSignalBeh',stopSignalBeh,'ttm',ttm,'ttx_history',ttx_history);
    
end
