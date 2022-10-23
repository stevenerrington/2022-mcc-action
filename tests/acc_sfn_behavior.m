% Code: Plot stopping behavior slowing
% 2022-10-16 | S P Errington

%% Analysis: Extract post-action slowing estimates.
% - For each behavioral session

ssd_cumul = [];
pnc_cumul = [];
monkey_cumul = [];
for session_i = 1:length(dataFiles_beh)
    
    % Setup workspace and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - clear loop variables
    clear rt_array
    
    % - get session information
    beh_idx = find(strcmp(dajo_datamap_curated.sessionBeh,dataFiles_beh(session_i)),1);
    sessionBeh = dajo_datamap_curated.sessionBeh(beh_idx);
    monkey = dajo_datamap_curated.monkey(beh_idx);
    
    rt_array = behavior(session_i).trialEventTimes.saccade-behavior(session_i).trialEventTimes.target;
    rt_array(rt_array < 100 | rt_array > 1000) = NaN;
    
    inh_function_session{session_i,1} = behavior(session_i).stopSignalBeh.inh_weibull.y;
    
    ssrt_session = behavior(session_i).stopSignalBeh.ssrt.integrationWeighted;
    
    stopBeh_plot(session_i,:) = table(sessionBeh,monkey,ssrt_session);
    
    RTdist.nostop{session_i,1} = cumulDist(rt_array(behavior(session_i).ttx.nostop.all.all));
    RTdist.noncanc{session_i,1} = cumulDist(rt_array(behavior(session_i).ttx.noncanceled.all.all));
    
    ssd_cumul = [ssd_cumul; behavior(session_i).stopSignalBeh.inh_SSD];
    pnc_cumul = [pnc_cumul; behavior(session_i).stopSignalBeh.inh_pnc'];
    monkey_cumul = [monkey_cumul; repmat(monkey,length(behavior(session_i).stopSignalBeh.inh_SSD),1)];
    
end

%% Figure: Inhibition function

inh_function_fig = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
clear inh_function_figure
inh_function_figure(1,1)= gramm('x',1:600,'y',inh_function_session);
inh_function_figure(1,1).stat_summary();
inh_function_figure(2,1)= gramm('x',ssd_cumul,'y',pnc_cumul);
inh_function_figure(1,1).stat_summary();
inh_function_figure(2,1).geom_point('alpha',0.2);

inh_function_figure.axe_property('YLim',[0 1]);
inh_function_figure(1,1).facet_grid([],stopBeh_plot.monkey);
inh_function_figure(2,1).facet_grid([],monkey_cumul);
inh_function_figure.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.figures,'sfn_inh_function.pdf');
set(inh_function_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(inh_function_fig,filename,'-dpdf') % then print it
close(inh_function_fig)

%% Figure: RT distributin

rt_dist_fig = figure('Renderer', 'painters', 'Position', [100 100 600 250]);
hold on
clear rt_dist_figure

for session_i = 1:length(dataFiles_beh)
    
    if strcmp(stopBeh_plot.monkey(session_i), 'dar')
        subplot(1,2,1); hold on
        
    else
        subplot(1,2,2); hold on
    end
    plot(RTdist.nostop{session_i,1}(:,1), RTdist.nostop{session_i,1}(:,2),'color',[colors.nostop 0.25])
    plot(RTdist.noncanc{session_i,1}(:,1), RTdist.noncanc{session_i,1}(:,2),'color',[colors.noncanc 0.25])
    xlim([100 700]); ylim([0 1])
end


% Once we're done with a page, save it and close it.
filename = fullfile(dirs.figures,'sfn_rt_dist.pdf');
set(rt_dist_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(rt_dist_fig,filename,'-dpdf') % then print it
close(rt_dist_fig)

%% Figure: SSRT histogram

ssrt_dist_fig = figure('Renderer', 'painters', 'Position', [100 100 300 250]);
hold on
clear ssrrt_dist_figure
histogram(stopBeh_plot.ssrt_session(strcmp(stopBeh_plot.monkey, 'dar')),60:5:180,'LineStyle','None')
hold on
histogram(stopBeh_plot.ssrt_session(strcmp(stopBeh_plot.monkey, 'jou')),60:5:180,'LineStyle','None')
xlim([60 180]); ylim([0 10])

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.figures,'sfn_ssrt_dist.pdf');
set(ssrt_dist_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(ssrt_dist_fig,filename,'-dpdf') % then print it
close(ssrt_dist_fig)

