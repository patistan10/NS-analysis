%% load different dataouts FOR NATURAL SCENES
clear all

i=1;
NSdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2209_2L_180302_suite2p_processed\processed_suite2p\2209_2L_180302_analysis\NatScenes\2209_2L_180302_NatScenes_dataOut.mat');i=i+1;
NSdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2210_NC_180228_suite2p_processed\processed_suite2p\2210_NC_180228_analysis\NatScenes\2210_NC_180228_NatScenes_dataOut.mat');i=i+1;
NSdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2212_NC_180227_suite2p_processed\processed_suite2p\2212_NC_180227_analysis\NatScenes\2212_NC_180227_NatScenes_dataOut.mat');i=i+1;
NSdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2231_1R1L_180304_suite2p_processed\processed_suite2p\2231_1R1L_180304_analysis\NatScenes\2231_1R1L_180304_NatScenes_dataOut.mat');i=i+1;

%get dimensions of all response matrices
dimensions = zeros(length(NSdataOuts),3);
for s=1:length(NSdataOuts)
    dimensions(s,1) = size(NSdataOuts{s}.dataOut.responseMatrix_1,1); %reps
    dimensions(s,2) = size(NSdataOuts{s}.dataOut.responseMatrix_1,2); %cells
    dimensions(s,3) = size(NSdataOuts{s}.dataOut.responseMatrix_1,3); %stim
end

%create a large dataOut 
dataOut.responseMatrix_1 = nan(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3));
dataOut.isRemovedBlock = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.hasLocomotion = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.totalNumStimuli = dimensions(1,3);
dataOut.totalNumCells = sum(dimensions(:,2));
dataOut.stimulus_ID = 'NatScenes';
dataOut.stats.global.response_average_pval_fdr = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.responsive_cells_p001_fdr_average_index = false(sum(dimensions(:,2)),1);
dataOut.stats.global.response_ACTUAL_avg_vals = nan(sum(dimensions(:,2)),dimensions(1,3));

%fill in large dataOut
count = 0;
count2 = 0;
FIRST_CELL = 1;
for s=1:length(NSdataOuts)
    %fill in response matrix
    current = NSdataOuts{s}.dataOut.responseMatrix_1;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.responseMatrix_1(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in removed block
    current = NSdataOuts{s}.dataOut.isRemovedBlock;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.isRemovedBlock(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in locomotion
    current = NSdataOuts{s}.dataOut.hasLocomotion;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.hasLocomotion(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;

    %fill in pvalue for median responses (using blcok average)
    current = NSdataOuts{s}.dataOut.stats.global.response_average_pval_fdr;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_pval_fdr(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;
    
    %fill in index of if cell is significant (pval<0.01) or not
    current = NSdataOuts{s}.dataOut.stats.global.response_average_pval_fdr;
    for c = 1:NSdataOuts{s}.dataOut.totalNumCells
        count = count+1;
        c_resps_pval = current(c,:);
        if min(c_resps_pval) <0.01
            dataOut.stats.global.responsive_cells_p001_fdr_average_index(count) = 1;
        end            
    end
    
    %fill in average value across trials after removing bad blocks
    current = NSdataOuts{s}.dataOut.responseMatrix_1;
    current2 = NSdataOuts{s}.dataOut.isRemovedBlock;
    for c = 1:NSdataOuts{s}.dataOut.totalNumCells
        count2=count2+1;
        c_resps = squeeze(current(:,c,:));
        c_bad = squeeze(current2(:,c,:));
        c_resps(c_bad==1) = nan;
        c_resps_avg = nanmean(c_resps,1);
        dataOut.stats.global.response_ACTUAL_avg_vals(count2,:) = c_resps_avg;        
    end 
    
    FIRST_CELL = FIRST_CELL + dimensions(s,2);
end
allcells = [1:sum(dimensions(:,2))]';
dataOut.stats.global.responsive_cells_p001_fdr_average = allcells(dataOut.stats.global.responsive_cells_p001_fdr_average_index);

save('dataOut_NatScenes_POOLED.mat','dataOut');


%% load different dataouts FOR GRATINGS
clear all

i=1;
GdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2209_2L_180302_suite2p_processed\processed_suite2p\2209_2L_180302_analysis\Gratings\2209_2L_180302_Grating_dataOut.mat');i=i+1;
GdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2210_NC_180228_suite2p_processed\processed_suite2p\2210_NC_180228_analysis\Gratings\2210_NC_180228_Grating_dataOut.mat');i=i+1;
GdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2212_NC_180227_suite2p_processed\processed_suite2p\2212_NC_180227_analysis\Gratings\2212_NC_180227_Grating_dataOut.mat');i=i+1;
GdataOuts{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2231_1R1L_180304_suite2p_processed\processed_suite2p\2231_1R1L_180304_analysis\Gratings\2231_1R1L_180304_Grating_dataOut.mat');i=i+1;

%get dimensions of all response matrices
dimensions = zeros(length(GdataOuts),3);
for s=1:length(GdataOuts)
    dimensions(s,1) = size(GdataOuts{s}.dataOut.responseMatrix_1,1); %reps
    dimensions(s,2) = size(GdataOuts{s}.dataOut.responseMatrix_1,2); %cells
    dimensions(s,3) = size(GdataOuts{s}.dataOut.responseMatrix_1,3); %stim
end

%create a large dataOut 
dataOut.responseMatrix_1 = nan(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3));
dataOut.isRemovedBlock = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.hasLocomotion = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.totalNumStimuli = dimensions(1,3);
dataOut.totalNumCells = sum(dimensions(:,2));
dataOut.stimulus_ID = 'NatScenes';
dataOut.stats.global.response_average_pval_fdr = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.responsive_cells_p001_fdr_average_index = false(sum(dimensions(:,2)),1);
dataOut.stats.global.response_ACTUAL_avg_vals = nan(sum(dimensions(:,2)),dimensions(1,3));

%fill in large dataOut
count = 0;
count2 = 0;
FIRST_CELL = 1;
for s=1:length(GdataOuts)
    %fill in response matrix
    current = GdataOuts{s}.dataOut.responseMatrix_1;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.responseMatrix_1(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in removed block
    current = GdataOuts{s}.dataOut.isRemovedBlock;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.isRemovedBlock(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in locomotion
    current = GdataOuts{s}.dataOut.hasLocomotion;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.hasLocomotion(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;

    %fill in pvalue for median responses (using blcok average)
    current = GdataOuts{s}.dataOut.stats.global.response_average_pval_fdr;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_pval_fdr(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;
    
    %fill in index of if cell is significant (pval<0.01) or not
    current = GdataOuts{s}.dataOut.stats.global.response_average_pval_fdr;
    for c = 1:GdataOuts{s}.dataOut.totalNumCells
        count = count+1;
        c_resps_pval = current(c,:);
        if min(c_resps_pval) <0.01
            dataOut.stats.global.responsive_cells_p001_fdr_average_index(count) = 1;
        end            
    end
    
    %fill in average value across trials after removing bad blocks
    current = GdataOuts{s}.dataOut.responseMatrix_1;
    current2 = GdataOuts{s}.dataOut.isRemovedBlock;
    for c = 1:GdataOuts{s}.dataOut.totalNumCells
        count2=count2+1;
        c_resps = squeeze(current(:,c,:));
        c_bad = squeeze(current2(:,c,:));
        c_resps(c_bad==1) = nan;
        c_resps_avg = nanmean(c_resps,1);
        dataOut.stats.global.response_ACTUAL_avg_vals(count2,:) = c_resps_avg;        
    end 
    
    FIRST_CELL = FIRST_CELL + dimensions(s,2);
end
allcells = [1:sum(dimensions(:,2))]';
dataOut.stats.global.responsive_cells_p001_fdr_average = allcells(dataOut.stats.global.responsive_cells_p001_fdr_average_index);

save('dataOut_Gratings_POOLED.mat','dataOut');

%% Pool gauss fit and BW for all cells
clear all

i=1;
gaussFits{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2209_2L_180302_suite2p_processed\processed_suite2p\2209_2L_180302_analysis\Gratings\Gauss_fit\gaussFit_results.mat');i=i+1;
gaussFits{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2210_NC_180228_suite2p_processed\processed_suite2p\2210_NC_180228_analysis\Gratings\Gauss_fit\gaussFit_results.mat');i=i+1;
gaussFits{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2212_NC_180227_suite2p_processed\processed_suite2p\2212_NC_180227_analysis\Gratings\Gauss_fit\gaussFit_results.mat');i=i+1;
gaussFits{i} = load('H:\ProcessedDataArchive\Pati\NatScene2\2231_1R1L_180304_suite2p_processed\processed_suite2p\2231_1R1L_180304_analysis\Gratings\Gauss_fit\gaussFit_results.mat');i=i+1;

%piece all variables together 
all_1minCV = [];
fit_oriBW = [];
fit_prefOri = [];
fit_rsquare = [];
pref_ori_estimate = [];
pref_sf_estimate = [];
gaussFitResults ={};
gaussgof ={};


for s = 1:length(gaussFits)
    all_1minCV = [all_1minCV gaussFits{s}.all_1minCV']; %rotate to be like other variables
    fit_oriBW = [fit_oriBW gaussFits{s}.fit_oriBW];
    fit_prefOri = [fit_prefOri gaussFits{s}.fit_prefOri];
    fit_rsquare = [fit_rsquare gaussFits{s}.fit_rsquare];
    pref_ori_estimate = [pref_ori_estimate gaussFits{s}.pref_ori_estimate];
    pref_sf_estimate = [pref_sf_estimate gaussFits{s}.pref_sf_estimate];
    gaussFitResults = [gaussFitResults gaussFits{s}.gaussFitResults];
    gaussgof = [gaussgof gaussFits{s}.gaussgof];
end

clear i s gaussFits
save('gaussFit_results_POOLED.mat');











