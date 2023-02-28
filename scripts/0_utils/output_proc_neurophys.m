function output_proc_neurophys(raster_in,filename_neural,saveDir)

names = fieldnames( raster_in );
subStr = 'DSP';
DSPstruct = rmfield( raster_in, names( find( cellfun( @isempty, strfind( names , subStr ) ) ) ) );
DSPnames = fieldnames(DSPstruct);

alignment_events = fieldnames(raster_in.(DSPnames{1}));

for DSPidx = 1:length(DSPnames)
    DSPlabel = DSPnames{DSPidx};
    
    for alignment_i = 1:length(alignment_events)
        raster = raster_in.(DSPlabel).(alignment_events{alignment_i});
        save(fullfile(saveDir,[filename_neural '_' DSPlabel '_' alignment_events{alignment_i} '_SPK']),'raster')
    end
    
end

end