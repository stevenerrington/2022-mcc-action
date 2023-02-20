% Code: Plot stopping behavior slowing
% 2022-10-16 | S P Errington

%% Analysis: Extract stopping behavior across sessions.
% - For each behavioral session

ssd_cumul = []; pnc_cumul = []; monkey_cumul = [];
rt_nostop_quantile = {}; rt_noncanc_quantile = {};
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
    
    rt_nostop_quantile{session_i,1} = quantile(RTdist.nostop{session_i,1}(:,1), [0.1:0.1:0.9]);
    rt_noncanc_quantile{session_i,1} = quantile(RTdist.noncanc{session_i,1}(:,1), [0.1:0.1:0.9]);
    
end


%% Figure: Generate figure

behavior_figure_out = figure('Renderer', 'painters', 'Position', [100 100 800 400]);
clear behavior_figure

% Inhibition function
behavior_figure(1,1)= gramm('x',1:600,'y',inh_function_session);
behavior_figure(1,1).geom_line('alpha',0.05);
behavior_figure(1,1).stat_summary();
behavior_figure(1,1).axe_property('XLim',[0 500],'YLim',[0 1]);
behavior_figure(1,1).facet_grid(stopBeh_plot.monkey,[]);
behavior_figure(1,1).no_legend();

% RT distribution
behavior_figure(1,2)= gramm('x',0.1:0.1:0.9,...
    'y',[rt_nostop_quantile;rt_noncanc_quantile],...
    'color',[repmat({'1_nostop'},length(rt_nostop_quantile),1);repmat({'2_noncanc'},length(rt_noncanc_quantile),1)]);
behavior_figure(1,2).stat_summary('geom',{'point','line','black_errorbar'});
behavior_figure(1,2).facet_grid(repmat(stopBeh_plot.monkey,2,1),[]);
behavior_figure(1,2).set_color_options('map',[colors.nostop; colors.noncanc]);
behavior_figure(1,2).axe_property('XLim',[0 1]);
behavior_figure(1,2).no_legend();

% SSRT histogram
behavior_figure(1,3)= gramm('x',stopBeh_plot.ssrt_session);
behavior_figure(1,3).stat_bin('edges',70:10:200,'fill','transparent');
behavior_figure(1,3).facet_grid(stopBeh_plot.monkey,[]);
behavior_figure(1,3).axe_property('XLim',[70 200]);
behavior_figure(1,3).axe_property('YLim',[0 15]);

behavior_figure.draw();



