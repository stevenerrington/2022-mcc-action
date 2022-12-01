

dajo_datamap_curated_DMFC = data_sessionCurate...
    (dajo_datamap,...
    'area', {'DMFC'}, 'monkey', {'dar','jou'}, 'signal', {'SPK'}, 'spacing', [50, 100, 150]);



% Get session information and indices
n_sessions = size(dajo_datamap_curated_DMFC,1);
dmfc_map_info = [];

dataFiles_beh = unique(dajo_datamap_curated_DMFC.sessionBeh);
dataFiles_neural = unique(dajo_datamap_curated_DMFC.dataFilename);

behavior = acc_stopping_extractBeh(dirs,dataFiles_beh);



for session_i = 1:n_sessions
    neuralFilename = dajo_datamap_curated_DMFC.dataFilename{session_i};
    behFilename = data_findBehFile(neuralFilename);
    beh_index = util_find_beh_index(behavior,behFilename);
    logInfo(session_i,:) = util_getLogInfo(neuralFilename);
    
    fprintf(['Extracting data for ' neuralFilename ': session %i of %i     \n'],...
        session_i, n_sessions);
        
    spk_data = load(fullfile(dirs.data,[neuralFilename '-spk.mat']));
    
    n_neurons= size(dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo,1);
    map_info = table();
    
    for neuron_i = 1:n_neurons
        site = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.site(neuron_i);
        session = dajo_datamap_curated_DMFC.dataFilename(session_i);
        monkey = dajo_datamap_curated_DMFC.monkey(session_i);
        unit = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.unitDSP(neuron_i);
        ap = logInfo.ap_stereo(session_i);
        ml = logInfo.ml_stereo(session_i);
        mua = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.flag_mua(neuron_i);
        noise = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.flag_noise(neuron_i);

        spk_width = util_getSpkWidth(spk_data,dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.unitWAV{neuron_i});
        
        map_info(neuron_i,:) = table(session,monkey,unit,site,ap,ml,mua,noise,spk_width);
    end
    
    dmfc_map_info = [dmfc_map_info; map_info];
end

% Clear up unknown areas and remove noise clusters
dmfc_map_info(dmfc_map_info.noise == 1,:) = [];


%%
spk_width_cutoff = 250;


%%
dmfc_session_list = unique(dmfc_map_info.session);

parfor (session_i = 1:length(dmfc_session_list),4)

    % We get the (neural) filename of the record of interest
    neuralFilename = dmfc_session_list{session_i};
    fprintf('Extracting: %s ... [%i of %i]  \n',neuralFilename,length(dmfc_session_list))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    import_data = struct(); import_data = load_spkFile(dirs,neuralFilename);
    
    % Convolve spike times to get continous trace
    spk_data_sdf = [];spk_data_spikes = [];
    [spk_data_sdf, spk_data_spikes] =...
        spk_alignTrials(behavior(behaviorIdx).trialEventTimes(:,[3,5,6,7,9,10]),...
        import_data.time, [-1000 2000]);

    % Then split this data into individual channels
    neuron_labels = {};
    neuron_labels = dmfc_map_info.unit(strcmp(dmfc_map_info.session,neuralFilename));
    
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


%% 
