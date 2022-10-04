function psd_analysis = lfp_getSessionPSD(alignedLFP, window)

% Clear and update workspace
clear CSDanalysis psd_analysis
fprintf(['Calculating CSD...\n'])

% Get alignment and channel information
alignName = 'reward';
channelNames = fieldnames(alignedLFP);

% For each channel in the penetration
for channelIdx = 1:length(channelNames)
    channel = channelNames{channelIdx};
    
    % Align and sort the LFP in the appropriate manner for the analysis
    CSDanalysis.(alignName).CSDarray(:,:,channelIdx) =...
        alignedLFP.(channel).(alignName)(:,:);
    CSDanalysis.(alignName).linearLFP(channelIdx,:) =...
        nanmean(CSDanalysis.(alignName).CSDarray(:,:,channelIdx));
end

% Once we have the data, reorient the array in the correct order for the
% analysis
CSDanalysis.(alignName).CSDarray =...
    permute(CSDanalysis.(alignName).CSDarray, [3 2 1]);
CSDanalysis.(alignName).CSDarray =...
    CSDanalysis.(alignName).CSDarray(:,window,:);

% Run the CSD suite (Westerberg, 2022)
CSDanalysis.(alignName).all =...
    SUITE_LAM(CSDanalysis.(alignName).CSDarray);

% Output PSD information
psd_analysis.psd = CSDanalysis.reward.all.PSD;
psd_analysis.psd_freq = CSDanalysis.reward.all.PSD_F;
psd_analysis.psd_norm = CSDanalysis.reward.all.PSD_NORM;


end



