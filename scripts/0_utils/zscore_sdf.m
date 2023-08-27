function z_sdf = zscore_sdf(sdf_in,baseline_win,event_alignment)


mean_fr = nanmean(nanmean(sdf_in.target(:,1000+baseline_win)));
std_fr = std(nanmean(sdf_in.target(:,1000+baseline_win)));

for trial_i = 1:size(sdf_in.target,1)
    z_sdf(trial_i,:) = (sdf_in.(event_alignment)(trial_i,:)-mean_fr)./std_fr;
end

end