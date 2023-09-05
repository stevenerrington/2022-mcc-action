function main_sdf_fig_out = main_sdf_figure(dirs,colors,...
     mcc_analysis_table,raw_neuron_idx,behavior,dataFiles_beh,...
     clusterNeurons, cluster_i, canc_sdf_in, nostop_sdf_in)


clear main_sdf_fig_out
main_sdf_fig_out = [];
main_sdf_fig_out = get_exampleNeuron_subplot(dirs,colors,...
     mcc_analysis_table,raw_neuron_idx,behavior,dataFiles_beh);

peak_fr_adj = [];
peak_fr_adj = max([max(canc_sdf_in(clusterNeurons{cluster_i},1000+[0:600]),[],2),...
   max(nostop_sdf_in(clusterNeurons{cluster_i},1000+[0:600]),[],2) ],[],2);
 
main_sdf_fig_out(3,1)=gramm('x',[-1000:2000],...
    'y',[num2cell(canc_sdf_in(clusterNeurons{cluster_i},:)./peak_fr_adj,2); num2cell(nostop_sdf_in(clusterNeurons{cluster_i},:)./peak_fr_adj,2)],...
    'color',[repmat({'1_canc'},length(clusterNeurons{cluster_i}),1);repmat({'2_nostop'},length(clusterNeurons{cluster_i}),1)]);
main_sdf_fig_out(3,1).stat_summary;
main_sdf_fig_out(3,1).geom_vline('xintercept',0);
main_sdf_fig_out(3,1).axe_property('XLim',[-250 750],'YLim',[0.5 1]);
main_sdf_fig_out(3,1).set_color_options('map',[colors.canceled;colors.nostop]);

main_sdf_fig_out(1,1).axe_property('XTick',[],'XColor',[1 1 1],'YTick',[],'YColor',[1 1 1]);
main_sdf_fig_out(2,1).axe_property('XTick',[],'XColor',[1 1 1]);
main_sdf_fig_out(3,1).axe_property('XTick',[-250 0 250 500 750]);


main_sdf_fig_out(1,1).set_layout_options('Position',[0.15 0.8 0.8 0.15],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.00 0.00],...
    'redraw',false);

main_sdf_fig_out(2,1).set_layout_options('Position',[0.15 0.45 0.8 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.00 0.00],...
    'redraw',false);

main_sdf_fig_out(3,1).set_layout_options('Position',[0.15 0.1 0.8 0.3],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.00 0.00],...
    'redraw',false);


end