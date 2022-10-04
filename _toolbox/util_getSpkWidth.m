function spk_width = util_getSpkWidth(spk_data,WAV_label)

wf_mean = nanmean(spk_data.spikes.waveform.(WAV_label));

width_nSamples = find(wf_mean == max(wf_mean)) - find(wf_mean == min(wf_mean));
fs = 24414.0625/1000;
spk_width = round(width_nSamples * fs);

end