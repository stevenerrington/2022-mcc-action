function dirs = get_dirs_mcc(user)

switch user
    case 'mac'
        dirs.root = '/Users/stevenerrington/Desktop/Projects/2022-mcc-action';
        dirs.toolbox = '/Users/stevenerrington/Desktop/Projects/toolbox';
        dirs.dajo_toolbox = '/Users/stevenerrington/Desktop/Projects/2022-dajo-toolbox';
        dirs.data = '/Volumes/Alpha/data/2021_Cmand_DaJo';
        
    case 'home'
        dirs.root = 'D:\projectCode\2022-mcc-action';
        dirs.toolbox = 'D:\projectCode\toolbox\';
        dirs.dajo_toolbox = 'D:\projectCode\2022-dajo-toolbox';
        dirs.data = 'D:\data\2021_Cmand_DaJo';
       
end

addpath(genpath(dirs.root));
addpath(genpath(dirs.toolbox));
addpath(genpath(dirs.dajo_toolbox));
addpath(genpath(dirs.data));

end

