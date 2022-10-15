acc_session_list = unique(acc_map_info.session);

parfor (session_i = 1:length(acc_session_list),4)

    % We get the (neural) filename of the record of interest
    neuralFilename = acc_session_list{session_i};
    fprintf('Extracting: %s ... [%i of %i]  \n',neuralFilename,length(acc_session_list))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    import_data = struct(); import_data = load_spkFile(dirs,neuralFilename);
    
    % Convolve spike times to get continous trace
    spk_data_sdf = [];spk_data_spikes = [];
    [spk_data_sdf, spk_data_spikes] =...
        spk_alignTrials(behavior(behaviorIdx).trialEventTimes(:,[3,5,6,7,9]),...
        import_data.time, [-1000 2000]);

    % Then split this data into individual channels
    neuron_labels = {};
    neuron_labels = acc_map_info.unit(strcmp(acc_map_info.session,neuralFilename));
    
    for neuron_i = 1:length(neuron_labels)
        SDF = []; Spikes = [];
        SDF = spk_data_sdf.(neuron_labels{neuron_i});
        Spikes = spk_data_spikes.(neuron_labels{neuron_i});
        
        save_filename_sdf = ...
            fullfile(dirs.root,'data','SDF',...
            [neuralFilename '_SDF_' neuron_labels{neuron_i} '.mat']);
        save_filename_spk = ...
            fullfile(dirs.root,'data','Spikes',...
            [neuralFilename '_Spikes_' neuron_labels{neuron_i} '.mat']);
        
        util_parsaveSDF(save_filename_sdf, SDF)
        util_parsaveSPK(save_filename_spk, Spikes)
    end
    
end
