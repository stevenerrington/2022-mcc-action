clear all; clc

% Set directories and workspace
root = 'C:\Users\Steven\Desktop\Projects\2022-acc-stopping\';
dataDir = 'S:\Users\Current Lab Members\Steven Errington\2021_DaJo\mat\';
cd(root)

% Load datamaps
load([root '\_data\2021-dajo-datamap.mat'])
ephysLog = importOnlineEphysLogMaster;

% Find ACC sessions
for sessionIdx = 1:size(dajo_datamap,1)
    acc_session_flag(sessionIdx,1) = sum(strcmp(dajo_datamap.neurophysInfo{sessionIdx}.area,'ACC')) > 0;
end
% Make an ACC specific dataframe to reference from
dajo_datamap_acc = dajo_datamap(acc_session_flag,:);

% Define ACC parameters
accParam.corticalDepth = 2000; % Cortical depth of ACC bank(in um);
accParam.sulcusDepth = 100; % Cortical depth of ACC bank(in um);

% Loop through ACC sessions
parfor sessionIdx = 1:size(dajo_datamap_acc,1)
    
    fprintf('Analysing session %i of %i | %s.          \n',...
        sessionIdx,size(dajo_datamap_acc,1),dajo_datamap_acc.session{sessionIdx});
    
    % Load in behavioural data
    behData = load([dataDir dajo_datamap_acc.behInfo(sessionIdx).dataFile]);
    
    % Get event times and trial types from the session for future alignments
    [ttx, ttx_history, trialEventTimes] =...
        processSessionTrials (behData.events.stateFlags_, behData.events.Infos_);
    
    % For each electrode in the session
    for elIdx = 1 : dajo_datamap_acc.nElectrodes(sessionIdx)
        
        % Load the Spike data and aligned it on SSRT (5) and Saccade (6)
        spkData = load([dataDir dajo_datamap_acc.neurophysInfo{sessionIdx}.spkFile{elIdx}]);
        tdtLFP = alignLFP(trialEventTimes(:,[5:6]), lfpData.lfp, [-500 1000]);
        
        % Find channels that cover dACC and vACC
        [accChannels] = findACCchannels(sessionIdx, elIdx, dajo_datamap_acc, ephysLog, accParam);
        
     
        
    end
end
