function [accChannels] = findACCchannels(sessionIdx, elIdx, dajo_datamap, ephysLog, accParam)

ephysLogIdx = find(strcmp(ephysLog.Session,...
    dajo_datamap.neurophysInfo{sessionIdx}.dataFilename{elIdx})); % Get ephysLog index for given penetration

acc_sulcus_ch = str2double(ephysLog.accBankChannel{ephysLogIdx}); % Get the channel which falls in ACC bank (visually inspected)

elSpacing = str2double(ephysLog.ElectrodeSpacing{ephysLogIdx}); % Get electrode spacing

acc_nChs = round(accParam.corticalDepth / elSpacing); % Find the number of channels covering the inputted expected ACC depth.
accSulcus_nChs = round(accParam.sulcusDepth / elSpacing); % Find the number of channels covering the inputted expected ACC depth.

accChannels.dACC = [acc_sulcus_ch - acc_nChs:1:acc_sulcus_ch-accSulcus_nChs]; % Find channels that occur within the dACC cortex size, up to the edge of the bank
accChannels.vACC = [acc_sulcus_ch + accSulcus_nChs:1:acc_sulcus_ch+acc_sulcus_ch]; % Find channels that occur within the vACC cortex size, up to the edge of the bank

accChannels.dACC(accChannels.dACC < 1 | accChannels.dACC > 32) = []; % Remove channels that fall out of range
accChannels.vACC(accChannels.vACC < 1 | accChannels.vACC > 32) = []; % Remove channels that fall out of range

end
