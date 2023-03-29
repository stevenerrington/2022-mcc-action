%{
Properties of midcingulate neurons in a saccade countermanding task
2022-10-11

Dependencies:
- gramm toolbox
%}




%% 1: Behavioral analysis
% Extract behavior across all relevant sessions (as identified in the
% curated dajo_datamap).
%   Behavior includes: stop behavior, RTs, value behavior.
behavior = mcc_stopping_extractBeh(dirs,dataFiles_beh);

sec1_figure_behaviorSummary

%% 2: 

