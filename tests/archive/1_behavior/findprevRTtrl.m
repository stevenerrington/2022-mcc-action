function prevrt_trls = findprevRTtrl(behavior)

% Define post-canceled no-stop trials (after first 10 trials - "warm up")
postc_trls = behavior.ttx_history.NS_after_C(behavior.ttx_history.NS_after_C > 10);

% Define all trials with a RT (non-canceled or no-stop)
rt_trls = sort([behavior.ttx.nostop.all.all;...
    behavior.ttx.noncanceled.all.all]);

% For each post-canceled no-stop trial, find the nearest previous trial
% with a RT
for trl = 1:length(postc_trls)
    [~, idx] = min(abs(rt_trls(rt_trls < postc_trls(trl)) -postc_trls(trl)));
    prevrt_trls(trl,1) = rt_trls(idx);
end

end
