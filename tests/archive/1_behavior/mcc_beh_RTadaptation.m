% Code: post-stopping slowing
% 2022-06-05 | S P Errington

%% Analysis: Extract post-action slowing estimates.
% - For each behavioral session
clear stopRTslowing
for session_i = 1:length(dataFiles_beh)
   
% Setup workspace and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% - clear loop variables 
clear RTarray trialHistory deltaRT

% - get session information
beh_idx = find(strcmp(dajo_datamap_curated.sessionBeh,dataFiles_beh(session_i)),1);
sessionBeh = dajo_datamap_curated.sessionBeh(beh_idx);
monkey = dajo_datamap_curated.monkey(beh_idx);

% - get session response latencies
RTarray = behavior(session_i).trialEventTimes.saccade - ...
    behavior(session_i).trialEventTimes.target;

% - find trials around the canceled trial with a response time
%    - note here, I've removed the first 10 trials, to allow the monkey to
%    settle and reduce variability
trialHistory.trl_after_c = behavior(session_i).ttx_history.NS_after_C...
    (behavior(session_i).ttx_history.NS_after_C > 10);
trialHistory.trl_before_c = findprevRTtrl(behavior(session_i));

% Analysis approaches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approach 1: get delta RT around the canceled trial
% - look at the change in RT in the last valid trial before the canceled
% trial, and the next valid trial with and RT. 
deltaRT = RTarray(trialHistory.trl_before_c) - RTarray(trialHistory.trl_after_c);
% - get the median change in RT. Negative values represent an RT slowing.
median_deltaRT = median(deltaRT);
% - determine whether this change is significant (Wilcoxon Sign)
[p_deltaRT, h_deltaRT] = signrank(deltaRT);

% Approach 2: compare differences between no-stop RT after canceled & no-stop trials
% Find the median RT on no-stop trials directly after canceled trials
median_postC = median(RTarray(behavior(session_i).ttx_history.NS_after_C));
% Find the median RT on no-stop trials directly after no-stop trials
median_postNS = median(RTarray(behavior(session_i).ttx_history.NS_after_NS));
median_postNC = median(RTarray(behavior(session_i).ttx_history.NS_after_NC));

% Output data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopRTslowing(session_i,:) = ... % Make a new row for this session
    table(sessionBeh, monkey,... % Include the key admin
    median_deltaRT, p_deltaRT, h_deltaRT,... % Info from approach 1
    median_postC, median_postNS,median_postNC); % Info from approach 2

end


%% Figure: Raw next trial RT, outcome-dependent

plotData = []; labelData = [];

plotData=...
    [stopRTslowing.median_postNS;...
    stopRTslowing.median_postNC;...
    stopRTslowing.median_postC];

labelData=...
    [repmat({'1_PostNS'},length(dataFiles_beh),1);...
    repmat({'2_PostNC'},length(dataFiles_beh),1);...
    repmat({'3_PostC'},length(dataFiles_beh),1)];

rt_slowing_fig = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
clear rt_slowing_figure
rt_slowing_figure(1,1)= gramm('x',labelData,'y',plotData,'color',labelData);
rt_slowing_figure(1,1).geom_jitter('alpha',0.2);
rt_slowing_figure(1,1).stat_summary('geom',{'point','line','black_errorbar'});
rt_slowing_figure.set_color_options('map',[colors.nostop;colors.noncanc;colors.canceled]);
rt_slowing_figure.axe_property('YLim',[250 450]);
rt_slowing_figure.facet_grid([],repmat(stopRTslowing.monkey,3,1));
rt_slowing_figure.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.figures,'rt_slowing_fig.pdf');
set(rt_slowing_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(rt_slowing_fig,filename,'-dpdf') % then print it
close(rt_slowing_fig)

%% Figure: Raw next trial RT, outcome-dependent

plotData = []; labelData = [];

plotData=...
    [stopRTslowing.median_deltaRT];

labelData=...
    [repmat({'1_PostC'},length(dataFiles_beh),1)];

alt_rt_slowing_fig = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
clear alt_rt_slowing_figure
alt_rt_slowing_figure(1,1)= gramm('x',labelData,'y',plotData,'color',labelData);
alt_rt_slowing_figure(1,1).geom_jitter('alpha',0.2);
alt_rt_slowing_figure(1,1).stat_summary('geom',{'point','line','black_errorbar'});
alt_rt_slowing_figure.set_color_options('map',[colors.canceled]);
alt_rt_slowing_figure.axe_property('YLim',[-50 25]);
alt_rt_slowing_figure.facet_grid([],stopRTslowing.monkey);
alt_rt_slowing_figure.draw();

% Once we're done with a page, save it and close it.
filename = fullfile(dirs.figures,'alt_rt_slowing_fig.pdf');
set(alt_rt_slowing_fig,'PaperSize',[20 10]); %set the paper size to what you want
print(alt_rt_slowing_fig,filename,'-dpdf') % then print it
close(alt_rt_slowing_fig)
