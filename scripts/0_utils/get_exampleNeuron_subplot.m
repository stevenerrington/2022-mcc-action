function example_neuron_subfig = get_exampleNeuron_subplot(dirs,colors,...
    mcc_analysis_table,raw_neuron_idx,behavior,dataFiles_beh)

spike_time_data = load(fullfile(dirs.root,'data','Spikes',...
    [mcc_analysis_table.session{raw_neuron_idx}...
    '_Spikes_' mcc_analysis_table.unit{raw_neuron_idx}]));

sdf_time_data = load(fullfile(dirs.root,'data','SDF',...
    [mcc_analysis_table.session{raw_neuron_idx}...
    '_SDF_' mcc_analysis_table.unit{raw_neuron_idx}]));

neuralFilename = mcc_analysis_table.session{raw_neuron_idx};
neuronLabel = mcc_analysis_table.unit{raw_neuron_idx};
fprintf('Extracting: %s - %s ...   \n',neuralFilename,neuronLabel)

%... and find the corresponding behavior file index
behFilename = data_findBehFile(neuralFilename);
behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));

canc_trial = []; nostop_trial = [];
canc_trial = behavior(behaviorIdx,:).ttx.canceled.all.all;
nostop_trial = behavior(behaviorIdx,:).ttx.nostop.all.all;

canc_raster = []; nostop_raster = [];
canc_raster = spike_time_data.Spikes.stopSignal_artifical(canc_trial);
nostop_raster = spike_time_data.Spikes.stopSignal_artifical(nostop_trial);

% Generate figure
clear example_neuron_subfig
example_neuron_subfig(1,1)=gramm('x',[canc_raster;nostop_raster],...
    'color',[repmat({'1_canc'},length(canc_raster),1);repmat({'2_nostop'},length(nostop_raster),1)]);
example_neuron_subfig(1,1).geom_raster;
example_neuron_subfig(1,1).axe_property('XLim',[-250 750]);
example_neuron_subfig(1,1).set_color_options('map',[colors.canceled;colors.nostop]);

example_neuron_subfig(2,1)=gramm('x',[-1000:2000],...
    'y',[sdf_time_data.SDF.stopSignal_artifical(canc_trial,:);sdf_time_data.SDF.stopSignal_artifical(nostop_trial,:)],...
    'color',[repmat({'1_canc'},length(canc_raster),1);repmat({'2_nostop'},length(nostop_raster),1)]);
example_neuron_subfig(2,1).stat_summary;
example_neuron_subfig(2,1).axe_property('XLim',[-250 750],'YLim',[5 25]);
example_neuron_subfig(2,1).geom_vline('xintercept',[0 behavior(behaviorIdx,:).stopSignalBeh.ssrt.integrationWeighted]);
example_neuron_subfig(2,1).set_color_options('map',[colors.canceled;colors.nostop]);

end