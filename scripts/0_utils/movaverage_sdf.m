function [window_sdf, window_time] = movaverage_sdf(sdf_in, window_size,window_step)

sdf_time = [-1000:2000];
sdf_index = 1:length(sdf_time);

window_end = 0; count = 0;

while window_end < max(sdf_index)
    window_start = (window_size/2)+(count*window_step);
    window_idx = window_start+[-window_size/2:+window_size/2]+1;
    window_end = max(window_idx);
    
    count = count+1;
    
    window_sdf(:,count) = nanmean(sdf_in(:,window_idx),2);
    window_time(:,count) = sdf_time(median(window_idx));
end
