
% Get session information and indices
n_sessions = size(dajo_datamap_curated,1);
all_map_info = [];

for session_i = 1:n_sessions
    neuralFilename = dajo_datamap_curated.dataFilename{session_i};
    behFilename = data_findBehFile(neuralFilename);
    beh_index = util_find_beh_index(behavior,behFilename);
    logInfo(session_i,:) = util_getLogInfo(neuralFilename);
    
    fprintf(['Extracting data for ' neuralFilename ': session %i of %i     \n'],...
        session_i, n_sessions);
    
    acs_ch_mapping = util_getACCchannels(logInfo,session_i);
    
    spk_data = load(fullfile(dirs.data,[neuralFilename '-spk.mat']));
    
    n_neurons= size(dajo_datamap_curated.spkInfo(session_i,1).unitInfo,1);
    map_info = table();
    
    for neuron_i = 1:n_neurons
        site = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.site(neuron_i);
        session = dajo_datamap_curated.dataFilename(session_i);
        monkey = dajo_datamap_curated.monkey(session_i);
        unit = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.unitDSP(neuron_i);
        depth = acs_ch_mapping.um_depth_acs(site);
        area = acs_ch_mapping.area(site);
        ap = logInfo.ap_stereo(session_i);
        ml = logInfo.ml_stereo(session_i);
        mua = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.flag_mua(neuron_i);
        noise = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.flag_noise(neuron_i);

        spk_width = util_getSpkWidth(spk_data,dajo_datamap_curated.spkInfo(session_i,1).unitInfo.unitWAV{neuron_i});
        
        map_info(neuron_i,:) = table(session,monkey,unit,site,depth,area,ap,ml,mua,noise,spk_width);
    end
    
    all_map_info = [all_map_info; map_info];
end


%%
all_map_info(strcmp(all_map_info.area,'?'),:) = [];
all_map_info(all_map_info.noise == 1,:) = [];
spk_width_cutoff = 250;

figure('Renderer', 'painters', 'Position', [100 100 250 500]);hold on
histogram(all_map_info.depth,-3000:150:3000,'LineStyle','None')
vline(0,'r-'); xlabel('Depth relative to ACS (\mum)'); ylabel('Neurons')
set(gca,'XDir','Reverse','view',[90 -90])

figure('Renderer', 'painters', 'Position', [100 100 600 500]);hold on;
% Narrow, typical spikes
subplot(2,2,1);
histogram(all_map_info.depth(0 < all_map_info.spk_width & all_map_info.spk_width <= spk_width_cutoff),-3000:150:3000,'LineStyle','None')
xlabel('Narrow Spikes'); vline(0,'r-')
set(gca,'XDir','Reverse','view',[90 -90])

% Broad, typical spikes
subplot(2,2,3);
histogram(all_map_info.depth(0 < all_map_info.spk_width & all_map_info.spk_width > spk_width_cutoff),-3000:150:3000,'LineStyle','None')
ylabel('Typical Spikes'); xlabel('Broad Spikes'); vline(0,'r-')
set(gca,'XDir','Reverse','view',[90 -90])

% Narrow, reverse spikes
subplot(2,2,2);
histogram(all_map_info.depth(0 > all_map_info.spk_width & abs(all_map_info.spk_width) <= spk_width_cutoff),-3000:150:3000,'LineStyle','None')
vline(0,'r-')
set(gca,'XDir','Reverse','view',[90 -90])

% Broad, reverse spikes
subplot(2,2,4);
histogram(all_map_info.depth(0 > all_map_info.spk_width & abs(all_map_info.spk_width) > spk_width_cutoff),-3000:150:3000,'LineStyle','None')
vline(0,'r-')
ylabel('Reverse Spikes'); vline(0,'r-')
set(gca,'XDir','Reverse','view',[90 -90])

%%
