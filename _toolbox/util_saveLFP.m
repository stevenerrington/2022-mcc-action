function util_saveLFP(lfp_aligned,neuralFilename,outputFolder)

channelNames = fieldnames(lfp_aligned);
nChannels = length(channelNames);

for channelIdx = 1:nChannels
    channel = channelNames{channelIdx};
    out_filename = [neuralFilename '-' channel '.mat'];
    LFP = lfp_aligned.(channel);
    save(fullfile(outputFolder,out_filename),'LFP')
end

end