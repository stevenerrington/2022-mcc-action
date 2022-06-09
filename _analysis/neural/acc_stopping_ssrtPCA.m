count = 0;

for electrode_i = 1:length(dataFiles_neural)
    nUnits = size(signal_average_spk{electrode_i}.individual.stopSignal_artifical,1);
    
    for unit_i = 1:nUnits
        count = count + 1;
        ns_ssd1(count,:) = signal_average_spk{electrode_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 1}(1,:);
        ns_ssd2(count,:) = signal_average_spk{electrode_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 1}(2,:);
        ns_ssd3(count,:) = signal_average_spk{electrode_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 1}(3,:);
        
        c_ssd1(count,:) = signal_average_spk{electrode_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 2}(1,:);
        c_ssd2(count,:) = signal_average_spk{electrode_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 2}(2,:);
        c_ssd3(count,:) = signal_average_spk{electrode_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 2}(3,:);
        
    end
    
end



%% PCA Analysis: Stop-Signal
% Set data parameters & windows
ephysWin = [-1000 2000]; winZero = abs(ephysWin(1));
plotWin = [-250:500]; 
analyseTime = [-100:200];
getColors

% Set PCA parameters
samplingRate = 1/1000; % inherent to the data. Do not change
numPCs = 8; % pick a number that will capture most of the variance
timeStep = 10; % smaller values will yield high resolution at the expense of computation time, default will sample at 20ms
withinConditionsOnly = false; % if true, will only analyze tanlging for times within the same condition

clear PCA_mainInput PCA_mainOutput
% Setup data for PCA analysis
PCA_mainInput(1).A = c_ssd1';
PCA_mainInput(1).times = plotWin';
PCA_mainInput(1).analyzeTimes = analyseTime';
PCA_mainInput(1).condition = 'SSD1 - canceled';
PCA_mainInput(1).color = [colors.canceled 0.25];

PCA_mainInput(2).A = c_ssd2';
PCA_mainInput(2).times = plotWin';
PCA_mainInput(2).analyzeTimes = analyseTime';
PCA_mainInput(2).condition = 'SSD2 - canceled';
PCA_mainInput(2).color = [colors.canceled 0.50];

PCA_mainInput(3).A = c_ssd3';
PCA_mainInput(3).times = plotWin';
PCA_mainInput(3).analyzeTimes = analyseTime';
PCA_mainInput(3).condition = 'SSD3 - canceled';
PCA_mainInput(3).color = [colors.canceled 1.00];

PCA_mainInput(4).A = ns_ssd1';
PCA_mainInput(4).times = plotWin';
PCA_mainInput(4).analyzeTimes = analyseTime';
PCA_mainInput(4).condition = 'SSD1 - canceled';
PCA_mainInput(4).color = [colors.nostop 0.25];

PCA_mainInput(5).A = ns_ssd2';
PCA_mainInput(5).times = plotWin';
PCA_mainInput(5).analyzeTimes = analyseTime';
PCA_mainInput(5).condition = 'SSD2 - canceled';
PCA_mainInput(5).color = [colors.nostop 0.50];

PCA_mainInput(6).A = ns_ssd3';
PCA_mainInput(6).times = plotWin';
PCA_mainInput(6).analyzeTimes = analyseTime';
PCA_mainInput(6).condition = 'SSD3 - canceled';
PCA_mainInput(6).color = [colors.nostop 1.00];
 

% Run PCA analysis
[~, PCA_mainOutput] = tangleAnalysis(PCA_mainInput, samplingRate,'numPCs',numPCs,'softenNorm',5 ,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % soft normalize neural data
pca_figure(PCA_mainInput,PCA_mainOutput)

