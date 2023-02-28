parfor (neural_i = 1:length(dataFiles_neural),4)
    
    % Get relevant file names
    filename_neural = dataFiles_neural{neural_i};
    filename_beh = data_findBehFile(filename_neural);
    
    % Update user
    fprintf('Extracting data from penetration %i of %i    | %s   \n',...
        neural_i, length(dataFiles_neural), filename_neural)
    
    % Load neural datafile
    data_in = load(fullfile(dirs.data, [filename_neural '-spk.mat']));
    
    % Define the alignment times for the session
    beh_idx = find(strcmp(behavior.sessionName,filename_beh));
    trialEventTimes = []; trialEventTimes = behavior.trialEventTimes{beh_idx};
    trialEventTimes = trialEventTimes(:,[3,5,6,9]);
    % 3: target, 5: saccade, 6: tone, 9: stop-signal
    
    % Loop through all neurons in a given session, and get the trial-by-trial sdf and raster plot  
    [raster] = spk_getRaster(trialEventTimes, data_in.spikes.time, [-1000 2000]);
    
    % Output all of these to a storage directory to be called later
     output_proc_neurophys(raster,filename_neural,...
         fullfile(dirs.root,'data','spk_rasters'))
end