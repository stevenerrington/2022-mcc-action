function data_in = util_clearFaultyCh(data_in, neuralFilename, ephysLog)

ephyslog_i = find(strcmp(neuralFilename,ephysLog.Session));
fault_ch = []; fault_ch = str2num(ephysLog.FaultCh{ephyslog_i});
nonfault_ch = []; nonfault_ch = find(~ismember(1:32,fault_ch));

if isempty(fault_ch)
    for fault_ch_i = 1:length(fault_ch)
        fault_ch_x = fault_ch(fault_ch_i);
        
        ch_label = ['LFP_' int2str(fault_ch_x)];
        data_in.data.(ch_label) = NaN(1,length(data_in.data.(ch_label)));
        
    end
end