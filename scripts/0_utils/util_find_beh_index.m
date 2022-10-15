function index = util_find_beh_index(behavior,behFilename)

for beh_i = 1:length(behavior)
    beh_file_list{beh_i,1} = behavior(beh_i).sessionName;
end

index = find(strcmp(beh_file_list,behFilename));

end