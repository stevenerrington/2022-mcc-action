
%% Extract: Example neuron data

example_session_file = 'jou-cmand1DR-ACC-20210621';
example_session_unit = 'DSP31b';
example_session_idx = find(strcmp(mcc_map_info.session,example_session_file) &...
    strcmp(mcc_map_info.unit,example_session_unit));

behFilename = data_findBehFile(example_session_file);
behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));

data_in_sdf = load(fullfile(dirs.root,'data','SDF',...
    [example_session_file '_SDF_' example_session_unit '.mat']));

data_in_spk = load(fullfile(dirs.root,'data','Spikes',...
    [example_session_file '_Spikes_' example_session_unit '.mat']));

sdf_noncanc_saccade = []; sdf_nostop_saccade = [];

trl_noncanc = []; trl_noncanc = behavior.ttx(behaviorIdx).noncanceled.all.all;
trl_nostop = []; trl_nostop = behavior.ttx(behaviorIdx).nostop.all.all;

example_sdf.saccade.noncanceled = num2cell(data_in_sdf.SDF.saccade(trl_noncanc,:),2);
example_sdf.saccade.nostop = num2cell(data_in_sdf.SDF.saccade(trl_nostop,:),2);
example_sdf.target.noncanceled = num2cell(data_in_sdf.SDF.target(trl_noncanc,:),2);
example_sdf.target.nostop = num2cell(data_in_sdf.SDF.target(trl_nostop,:),2);

spk_times_noncanc_target = {}; spk_times_nostop_target = {};
spk_times_noncanc_saccade = {}; spk_times_nostop_saccade = {};

for trl_i = 1:length(trl_noncanc)
    spk_times_noncanc_target{trl_i,1} = data_in_spk.Spikes.target{trl_noncanc(trl_i)};
    spk_times_noncanc_saccade{trl_i,1} = data_in_spk.Spikes.saccade{trl_noncanc(trl_i)};
end
for trl_i = 1:length(trl_nostop)
    spk_times_nostop_target{trl_i,1} = data_in_spk.Spikes.target{trl_nostop(trl_i)};
    spk_times_nostop_saccade{trl_i,1} = data_in_spk.Spikes.saccade{trl_nostop(trl_i)};
end



%% Example unit
clear error_pop_fig1

% Input data ---------------------------
error_pop_fig1(1,1) = gramm('x',[spk_times_nostop_target;spk_times_noncanc_target],...
    'color',[repmat({'1_nostop'},length(trl_nostop),1);...
    repmat({'2_noncanc'},length(trl_noncanc),1)]);

error_pop_fig1(2,1) = gramm('x',[spk_times_nostop_saccade;spk_times_noncanc_saccade],...
    'color',[repmat({'1_nostop'},length(trl_nostop),1);...
    repmat({'2_noncanc'},length(trl_noncanc),1)]);

error_pop_fig1(3,1)= gramm('x',-1000:2000,...
    'y',[example_sdf.target.nostop;example_sdf.target.noncanceled],...
    'color',[repmat({'1_nostop'},length(trl_nostop),1);...
    repmat({'2_noncanc'},length(trl_noncanc),1)]);

error_pop_fig1(4,1)= gramm('x',-1000:2000,...
    'y',[example_sdf.saccade.nostop;example_sdf.saccade.noncanceled],...
    'color',[repmat({'1_nostop'},length(trl_nostop),1);...
    repmat({'2_noncanc'},length(trl_noncanc),1)]);


%% Population SDF
% Input data ---------------------------
cluster_id = 4;
error_pop_fig1(5,1)= gramm('x',-1000:2000,...
    'y',[num2cell(norm_sdf.target.nostop(clusterNeurons{cluster_id},:),2);...
    num2cell(norm_sdf.target.noncanc(clusterNeurons{cluster_id},:),2)],...
    'color',[repmat({'1_nostop'},length(clusterNeurons{cluster_id}),1);...
    repmat({'2_noncanc'},length(clusterNeurons{cluster_id}),1)]);
error_pop_fig1(6,1)= gramm('x',-1000:2000,...
    'y',[num2cell(norm_sdf.saccade.nostop(clusterNeurons{cluster_id},:),2);...
    num2cell(norm_sdf.saccade.noncanc(clusterNeurons{cluster_id},:),2)],...
    'color',[repmat({'1_nostop'},length(clusterNeurons{cluster_id}),1);...
    repmat({'2_noncanc'},length(clusterNeurons{cluster_id}),1)]);

%% 
% Figure type


error_pop_fig1(1,1).geom_raster(); error_pop_fig1(2,1).geom_raster(); 
error_pop_fig1(3,1).stat_summary(); error_pop_fig1(4,1).stat_summary(); 
error_pop_fig1(5,1).stat_summary(); error_pop_fig1(6,1).stat_summary(); 

% Axe properties

error_pop_fig1(1,1).axe_property('XLim',[-400 200],'XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);
error_pop_fig1(2,1).axe_property('XLim',[-200 1800],'XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);


error_pop_fig1(3,1).axe_property('XLim',[-400 200],'YLim',[10 25],'XTick',[],'XColor',[1 1 1]);
error_pop_fig1(4,1).axe_property('XLim',[-200 1800],'YLim',[10 25],'XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);

error_pop_fig1(5,1).axe_property('XLim',[-400 200],'YLim',[-2 10]);
error_pop_fig1(6,1).axe_property('XLim',[-200 1800],'YLim',[-2 10],'YTick',[],'YColor',[1 1 1]);

% Colormaps
error_pop_fig1(1,1).set_color_options('map',[colors.nostop;colors.noncanc]);
error_pop_fig1(2,1).set_color_options('map',[colors.nostop;colors.noncanc]);
error_pop_fig1(3,1).set_color_options('map',[colors.nostop;colors.noncanc]);
error_pop_fig1(4,1).set_color_options('map',[colors.nostop;colors.noncanc]);
error_pop_fig1(5,1).set_color_options('map',[colors.nostop;colors.noncanc]);
error_pop_fig1(6,1).set_color_options('map',[colors.nostop;colors.noncanc]);

% Markers
error_pop_fig1(2,1).geom_vline('xintercept',[0 500 1000],'style','k-');
error_pop_fig1(4,1).geom_vline('xintercept',[0 500 1000],'style','k-');
error_pop_fig1(6,1).geom_vline('xintercept',[0 500 1000],'style','k-');



% Position figure: [left bottom width height]

error_pop_fig1(1,1).set_layout_options...
    ('Position',[0.1 0.8 0.15 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... 
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);
error_pop_fig1(2,1).set_layout_options...
    ('Position',[0.3 0.8 0.15*3.333 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... 
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

error_pop_fig1(3,1).set_layout_options...
    ('Position',[0.1 0.45 0.15 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... 
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);
error_pop_fig1(4,1).set_layout_options...
    ('Position',[0.3 0.45 0.15*3.333 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... 
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


error_pop_fig1(5,1).set_layout_options...
    ('Position',[0.1 0.1 0.15 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... 
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);
error_pop_fig1(6,1).set_layout_options...
    ('Position',[0.3 0.1 0.15*3.333 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... 
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

error_pop_fig1_out = figure('Renderer', 'painters', 'Position', [100 100 500 600]);
error_pop_fig1.draw
