function logInfo = util_getLogInfo(neuralFilename)
    ephysLog = importOnlineEphysLogMaster;
    log_i = find(strcmp(ephysLog.Session,neuralFilename));
    
    ap_stereo = str2num(ephysLog.AP_Stereotaxic{log_i});
    ml_stereo = str2num(ephysLog.ML_Stereotaxic{log_i});
    electrode_depth = str2num(ephysLog.ElectrodeSettleDepth{log_i});
    acs_ch = str2num(ephysLog.accBankChannel{log_i});
    electrode_spc = str2num(ephysLog.ElectrodeSpacing{log_i});
    
    if isempty(acs_ch)
        acs_ch = NaN;
    end
    
    neuralFilename = {neuralFilename};
    logInfo = table(neuralFilename,log_i,ap_stereo,ml_stereo,...
        electrode_depth,acs_ch,electrode_spc);
end
