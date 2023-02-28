function [raster] = spk_getRaster(trialEventTimes, spkTimes, timeWin)

names = fieldnames( spkTimes );
subStr = 'DSP';
DSPstruct = rmfield( spkTimes, names( find( cellfun( @isempty, strfind( names , subStr ) ) ) ) );
DSPnames = fieldnames(DSPstruct);

for DSPidx = 1:length(DSPnames)
    DSPlabel = DSPnames{DSPidx};
    
    eventNames = fieldnames(trialEventTimes);
    eventNames = eventNames(1:length(eventNames)-3);
    
    alignedRaster = {};
    
    for alignIdx = 1:length(eventNames)
        alignTimes = round(trialEventTimes.(eventNames{alignIdx})(:));
        
        alignedSDF_event = nan(length(alignTimes),range(timeWin)+1);
        raster_temp = zeros(length(alignTimes),range(timeWin)+1);
        
        spk_aligntemp = {};
        for trl_i = 1:length(alignTimes)
            if isnan(alignTimes(trl_i)) | alignTimes(trl_i) == 0
                spk_aligntemp{trl_i,1} = [];
            else
                spk_aligntemp{trl_i,1} = intersect(spkTimes.(DSPlabel),alignTimes(trl_i)+[timeWin(1):timeWin(2)])-...
                    alignTimes(trl_i);
                
                trl_spkTimes = spk_aligntemp{trl_i,1}+abs(min(timeWin(1)))+1;
                raster_temp(trl_i,trl_spkTimes) = 1;
            end
        end
        
        alignedRaster{alignIdx} = raster_temp;
        
    end
    
    for alignIdx = 1:length(eventNames)
        raster.(DSPlabel).(eventNames{alignIdx}) = alignedRaster{alignIdx};
    end
    
    clear alignedSDF aligned_spkTimes
end

end

