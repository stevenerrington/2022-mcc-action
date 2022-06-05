function signal_average = acc_stopping_extractSDF(dirs,dataFiles_beh,dataFiles_neural,behavior)

parfor dataFileIdx = 1:length(dataFiles_neural)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_neural{dataFileIdx},dataFileIdx,length(dataFiles_neural))
    
    % We then get the (neural) filename of the record of interest
    neuralFilename = dataFiles_neural{dataFileIdx};
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in data-file
    import_data = struct(); import_data = load_spkFile(dirs,neuralFilename);
    
    % Convolve spike times to get continous trace
    spk_data = [];
    spk_data = spk_alignTrials(behavior(behaviorIdx).trialEventTimes(:,[3,6,9,10]),...
        import_data.time, [-1000 2000]);

    ssd_in = [behavior(behaviorIdx).stopSignalBeh.midSSDidx - 1,...
        behavior(behaviorIdx).stopSignalBeh.midSSDidx,...
        behavior(behaviorIdx).stopSignalBeh.midSSDidx + 1];
    
    signal_average{dataFileIdx} = neural_ttmSignalAverage(spk_data,...
        behavior(behaviorIdx).ttm.C,...
        'units',fieldnames(spk_data),...
        'events',{'target','stopSignal_artifical','ssrt'},...
        'ssd',ssd_in,...
        'weighting',{'on'})
    
end