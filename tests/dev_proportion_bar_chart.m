
n_dMFC = sum(strcmp(dmfc_map_info.area, 'DMFC'));
n_dMCC = sum(strcmp(mcc_map_info.area, 'dMCC'));
n_vMCC = sum(strcmp(mcc_map_info.area, 'vMCC') | strcmp(mcc_map_info.area, 'ACS'));


% Action encoding
n_canc_trialType_dMFC = (sum(dmfc_analysis_table.glm_trial == 1 &...
    strcmp(dmfc_map_info.area, 'DMFC'))./n_dMFC)*100;

n_canc_trialType_dMCC = (sum(mcc_analysis_table.glm_trial == 1 &...
    strcmp(mcc_map_info.area, 'dMCC'))./n_dMCC)*100;

n_canc_trialType_vMCC = (sum(mcc_analysis_table.glm_trial == 1 &...
    strcmp(mcc_map_info.area, 'vMCC'))./n_vMCC)*100;

% SSD encoding
n_canc_ssd_dMFC = (sum(dmfc_analysis_table.glm_ssd_canc == 1 &...
    strcmp(dmfc_map_info.area, 'DMFC'))./n_dMFC)*100;

n_canc_ssd_dMCC = (sum(mcc_analysis_table.glm_ssd_canc == 1 &...
    strcmp(mcc_map_info.area, 'dMCC'))./n_dMCC)*100;

n_canc_ssd_vMCC = (sum(mcc_analysis_table.glm_ssd_canc == 1 &...
    strcmp(mcc_map_info.area, 'vMCC'))./n_vMCC)*100;

% Value encoding
n_canc_value_dMFC = (sum(dmfc_analysis_table.glm_value_canc == 1 &...
    strcmp(dmfc_map_info.area, 'DMFC'))./n_dMFC)*100;

n_canc_value_dMCC = (sum(mcc_analysis_table.glm_value_canc == 1 &...
    strcmp(mcc_map_info.area, 'dMCC'))./n_dMCC)*100;

n_canc_value_vMCC = (sum(mcc_analysis_table.glm_value_canc == 1 &...
    strcmp(mcc_map_info.area, 'vMCC'))./n_vMCC)*100;


%% Figure
clear figure_plot
figure_plot(1,1)=gramm('x',{'1_Trial Type','1_Trial Type','1_Trial Type','2_SSD','2_SSD','2_SSD','3_Value','3_Value','3_Value'},...
'y',[n_canc_trialType_dMFC,n_canc_trialType_dMCC, n_canc_trialType_vMCC, n_canc_ssd_dMFC, n_canc_ssd_dMCC, n_canc_ssd_vMCC, n_canc_value_dMFC, n_canc_value_dMCC, n_canc_value_vMCC],...
'color',{'1_dMFC','2_dMCC','3_vMCC','1_dMFC','2_dMCC','3_vMCC','1_dMFC','2_dMCC','3_vMCC'});
figure_plot(1,1).stat_summary('geom',{'bar'});
figure_plot(1,1).axe_property('YLim',[0 100]);
figure_plot(1,1).set_names('y','% significant neurons in area');
figure
figure_plot.draw