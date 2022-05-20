clear all; clc

% Set directories and workspace
root = 'C:\Users\Steven\Desktop\Projects\2022-acc-stopping\';
dataDir = 'S:\Users\Current Lab Members\Steven Errington\2021_DaJo\mat\';
cd(root)

% Load datamap
load([root '\_data\2021-dajo-datamap.mat'])

% Loop session
for session = 1:size(dajo_datamap,1)
    
    % Let the user know where we are at!
    fprintf('Extracting event time from %s | %i of %i \n',...
    dajo_datamap.behInfo(session,1).dataFile,...
    session, size(dajo_datamap,1))

    % Load beh data
    beh_data = load(fullfile(dataDir, dajo_datamap.behInfo(session,1).dataFile));
    
    % Get trial indices and key event times for alignment
    [ttx, ~, trialEventTimes] = processSessionTrials...
        (beh_data.events.stateFlags_,...
         beh_data.events.Infos_);
     
     % Get an estimate of stop-signal reaction time
     [stopSignalBeh, ~] = extractStopBeh...
         (beh_data.events.stateFlags_,...
         beh_data.events.Infos_,...
         ttx);
     
     % Find trials with no stop-signals
     no_stopSignal_trls = find(isnan(trialEventTimes.stopSignal));
     
     % For each no stop-signal trial
     for trlIdx = 1:length(no_stopSignal_trls)
         trl = no_stopSignal_trls(trlIdx);
         
         % Get an estimated SSRT time for trials where there was a target.
         trialEventTimes.ssrt(trl) = ...
         	trialEventTimes.target(trl) + ... % Get the target time (as it's NaN for no target, we won't get a value).
            round(beh_data.events.stateFlags_.LastSsdIdx(trl)*(1000/60)) +... % ... add the SSD (ms) from the previous stop trial
            stopSignalBeh.ssrt.integrationWeighted; % ... and add SSRT. 

     end
     
     writetable(trialEventTimes,...
         [root '\_data\event-times\' dajo_datamap.behInfo(session,1).dataFile '-events.csv']);

end

