




%% Spiking activity
signal_average_spk = acc_stopping_extractSDF(dirs,dataFiles_beh,dataFiles_neural,behavior);

signal_collapse = neural_collapseSignalSession(signal_average_spk,...
    'events',{'target','stopSignal_artifical','ssrt'},...
    'conditions',{'C','GO'},...
    'conditions_map',[1 2]);

acc_stopping_ssrtClustering
acc_stopping_ssrtPCA

