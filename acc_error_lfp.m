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
        
        % Load the LFP data and aligned it on SSRT (5) and Saccade (6)
        lfpData = load([dataDir dajo_datamap_acc.neurophysInfo{sessionIdx}.lfpFile{elIdx}]);
        tdtLFP = alignLFP(trialEventTimes(:,[5:6]), lfpData.lfp, [-500 1000]);
        
        % Find channels that cover dACC and vACC
        [accChannels] = findACCchannels(sessionIdx, elIdx, dajo_datamap_acc, ephysLog, accParam);
        
        % Get averaged event aligned LFPs for SSRT and Saccade
        % dACC
        if ~isnan(accChannels.dACC) & ~isempty(accChannels.dACC)
            for dACC_ch_idx = 1:length(accChannels.dACC)
                accLFP{sessionIdx}.dACC.error.noncanceled(dACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.dACC(dACC_ch_idx))]).saccade...
                    (ttx.noncanceled.all.all,:));
                accLFP{sessionIdx}.dACC.error.nostop(dACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.dACC(dACC_ch_idx))]).saccade...
                    (ttx.nostop.all.all,:));
                accLFP{sessionIdx}.dACC.ssrt.canceled(dACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.dACC(dACC_ch_idx))]).ssrt...
                    (ttx.canceled.all.all,:));
                accLFP{sessionIdx}.dACC.ssrt.nostop(dACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.dACC(dACC_ch_idx))]).ssrt...
                    (ttx.nostop.all.all,:));
            end
        end
        
        % vACC
        if ~isnan(accChannels.vACC) & ~isempty(accChannels.vACC)
            for vACC_ch_idx = 1:length(accChannels.vACC)
                accLFP{sessionIdx}.dACC.error.noncanceled(vACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.vACC(vACC_ch_idx))]).saccade...
                    (ttx.noncanceled.all.all,:));
                accLFP{sessionIdx}.dACC.error.nostop(vACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.vACC(vACC_ch_idx))]).saccade...
                    (ttx.nostop.all.all,:));
                accLFP{sessionIdx}.vACC.ssrt.canceled(vACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.vACC(vACC_ch_idx))]).ssrt...
                    (ttx.canceled.all.all,:));
                accLFP{sessionIdx}.vACC.ssrt.nostop(vACC_ch_idx,:) = ...
                    nanmean(tdtLFP.aligned.(['LFP_' int2str(accChannels.vACC(vACC_ch_idx))]).ssrt...
                    (ttx.nostop.all.all,:));
            end
        end
    end
end
