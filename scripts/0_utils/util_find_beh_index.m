function index = util_find_beh_index(behavior,behFilename)

for beh_i = 1:size(behavior,1)
    beh_file_list{beh_i,1} = behavior.sessionName(beh_i);
end

index = find(strcmp(beh_file_list,behFilename));

end