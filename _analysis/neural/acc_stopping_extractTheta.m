lfpDir = 'S:\Users\Current Lab Members\Steven Errington\analysis\2021-dajo-lfp\';
outDir = 'C:\Users\Steven\Desktop\Projects\2022-acc-stopping\_data\lfp\sessionTheta';

%% Main script
for dataFileIdx = 1:length(dataFiles_neural)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_neural{dataFileIdx},dataFileIdx,length(dataFiles_neural))
    
    % We then get the (neural) filename of the record of interest
    neuralFilename = dataFiles_neural{dataFileIdx};
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));

    thetaLFP = {};
    parfor channel_i = 1:32
        channel = ['LFP_' int2str(channel_i)];
        in_filename = [neuralFilename '-' channel '.mat'];
        data_in = load(fullfile(lfpDir,in_filename));
        
        thetaLFP{channel_i} = util_bandpassFilter(double(data_in.LFP.stopSignal_artifical),...
            1000, 1, filter.band, filter.label, 2);
        
    end
    
    outFile = fullfile(outDir, [neuralFilename '-stopSignal-theta']);
    save(outFile,'thetaLFP', '-v7.3')
    
end


%%
parfor (dataFileIdx = 1:length(dataFiles_neural),10)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_neural{dataFileIdx},dataFileIdx,length(dataFiles_neural))
    
    % We then get the (neural) filename of the record of interest
    neuralFilename = dataFiles_neural{dataFileIdx};
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    inFile = fullfile(outDir, [neuralFilename '-stopSignal-theta']);
    data_in = load(inFile);
    
    ssd_x = [behavior(behaviorIdx).stopSignalBeh.midSSDidx - 1,...
        behavior(behaviorIdx).stopSignalBeh.midSSDidx,...
        behavior(behaviorIdx).stopSignalBeh.midSSDidx + 1];
    
    for ssd_i = 1:3
        trl_in_c = []; trl_in_ns = [];  trl_in_nc = [];
        trl_in_c = behavior(behaviorIdx).ttm.C.C{ssd_x(ssd_i)};
        trl_in_ns = behavior(behaviorIdx).ttm.C.GO{ssd_x(ssd_i)};
        trl_in_nc = behavior(behaviorIdx).ttm.NC.NC{ssd_x(ssd_i)};
        
        for channel_i = 1:32
            tta_pwr_c{dataFileIdx,ssd_i}(channel_i,:) = nanmean(data_in.thetaLFP{channel_i}.theta_pwr(trl_in_c,:));
            tta_pwr_ns{dataFileIdx,ssd_i}(channel_i,:) = nanmean(data_in.thetaLFP{channel_i}.theta_pwr(trl_in_ns,:));
            tta_pwr_nc{dataFileIdx,ssd_i}(channel_i,:) = nanmean(data_in.thetaLFP{channel_i}.theta_pwr(trl_in_nc,:));
        end
        
    end
end


%%
for dataFileIdx = 1:length(dataFiles_neural)
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    ssrt = round(behavior(behaviorIdx).stopSignalBeh.ssrt.integrationWeighted);
    
    ssd_i = behavior(behaviorIdx).stopSignalBeh.midSSDidx;
    ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(ssd_i);
    
    
    thetaWindow_bl = []; thetaWindow_bl = 1000-ssd-200:1000-ssd;
    thetaWindow_stop = []; thetaWindow_stop = 1000:1000+ssrt;
    thetaWindow_ssrt = []; thetaWindow_ssrt = 1000+ssrt:1000+ssrt+ssrt;
    
 
    tta_pwr_mean_ns_bl(dataFileIdx,1) = nanmean(nanmean(tta_pwr_ns{dataFileIdx,2}(:,thetaWindow_bl)));
    tta_pwr_mean_nc_bl(dataFileIdx,1) = nanmean(nanmean(tta_pwr_nc{dataFileIdx,2}(:,thetaWindow_bl)));
    tta_pwr_mean_c_bl(dataFileIdx,1) = nanmean(nanmean(tta_pwr_c{dataFileIdx,2}(:,thetaWindow_bl)));

    
    tta_pwr_mean_ns_stopping(dataFileIdx,1) = nanmean(nanmean(tta_pwr_ns{dataFileIdx,2}(:,thetaWindow_stop)));
    tta_pwr_mean_nc_stopping(dataFileIdx,1) = nanmean(nanmean(tta_pwr_nc{dataFileIdx,2}(:,thetaWindow_stop)));
    tta_pwr_mean_c_stopping(dataFileIdx,1) = nanmean(nanmean(tta_pwr_c{dataFileIdx,2}(:,thetaWindow_stop)));

    tta_pwr_mean_ns_ssrt(dataFileIdx,1) = nanmean(nanmean(tta_pwr_ns{dataFileIdx,2}(:,thetaWindow_ssrt)));
    tta_pwr_mean_nc_ssrt(dataFileIdx,1) = nanmean(nanmean(tta_pwr_nc{dataFileIdx,2}(:,thetaWindow_ssrt)));
    tta_pwr_mean_c_ssrt(dataFileIdx,1) = nanmean(nanmean(tta_pwr_c{dataFileIdx,2}(:,thetaWindow_ssrt)));
end

data = [(tta_pwr_mean_ns_stopping./tta_pwr_mean_ns_bl)*100-100;...
    (tta_pwr_mean_nc_stopping./tta_pwr_mean_nc_bl)*100-100;...
    (tta_pwr_mean_c_stopping./tta_pwr_mean_c_bl)*100-100;...
    (tta_pwr_mean_ns_ssrt./tta_pwr_mean_ns_bl)*100-100;...
    (tta_pwr_mean_nc_ssrt./tta_pwr_mean_nc_bl)*100-100;...
    (tta_pwr_mean_c_ssrt./tta_pwr_mean_c_bl)*100-100];

monkey_labels = repmat(dajo_datamap_curated.monkey,6,1);

trial_labels = repmat([repmat({'1_no-stop'},length(dataFiles_neural),1);...
    repmat({'2_non-canc'},length(dataFiles_neural),1);...
    repmat({'3_canc'},length(dataFiles_neural),1)],2,1);

epoch_labels = [repmat({'1_stopping'},length(dataFiles_neural)*3,1);...
    repmat({'2_ssrt'},length(dataFiles_neural)*3,1)];
     

%% Figure
clear test
figure('Renderer', 'painters', 'Position', [100 100 300 300])
test(1,1)= gramm('x',trial_labels,...
    'y',data,'color',trial_labels,'subset',strcmp(epoch_labels,'1_stopping'));
test(1,1).stat_summary('geom',{'point','errorbar','lines'},'type','sem');
test.set_color_options('map',[colors.nostop;colors.noncanc;colors.canceled]);
test.axe_property('YLim',[15 35])
test.no_legend
test.draw();

% data = data*10e14;
%% ANOVA table
outtable = table(monkey_labels,trial_labels,epoch_labels,data);
writetable(outtable,fullfile('C:\Users\Steven\Desktop\Projects\2022-acc-stopping\_data\tables','jasp_thetapower_stopping.csv'),'WriteRowNames',true)