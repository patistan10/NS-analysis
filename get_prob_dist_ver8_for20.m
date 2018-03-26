function [RespMatrix,RespMatrixTraining,RespMatrixTesting,RespMatrixAllTrials,ProbDist,BinThresh,removedBlocks,selected_cells] = get_prob_dist_ver8_for20(ExpID,selected_cells,bins,folds,test_set)
%UNTITLED2 Summary of this function goes here
%   SINGLE SESSION decoding


    if isempty(ExpID)
        error('Please enter experiment ID')
    end

    %ExpID = '2162_1R2L_170906'; %%%%%%will be specified as input (delete after testing)
    current_dir = pwd;
    cd ../..
    load(sprintf('%s_NatScenes_dataOut.mat',ExpID));
    cd(current_dir)
    
    %dataOut.totalNumStimuli = 10; % FOR CONTRAST TEST
    %response_matrix_all = dataOut.responseMatrix_1; %response matrix for all trials, all cells, all stim
    %response_matrix_all = response_matrix_all(1:15,:,:); %used to test what happens if less data
    response_matrix_all = dataOut.responseMatrix_1(:,:,:); % FOR MIXTURE
    
    
    %removedBlocks = squeeze(dataOut.isRemovedBlock(:,1,:));
    %removedBlocks = squeeze(dataOut.isRemovedBlock(1:15,1,:));%used to test what happens if less data
    removedBlocks = squeeze(dataOut.isRemovedBlock(:,1,:)); % FOR MIXTURE
    
    
    if isempty(selected_cells)
         %selected_cells = (1:dataOut.totalNumCells);         
        selected_cells = dataOut.stats.global.responsive_cells_p005_fdr_average'; %cells responsive to at least one nat stim
        %selected_cells = selected_cells(91:120);
    end
    if isempty(bins)
        bins = 3; %default number of bins
    end
    if isempty(folds)
        folds = 10;
    end
    if isempty(test_set)
        test_set = 1;
    end
    
    RespMatrix = response_matrix_all(:,selected_cells,:); 
    ProbDist = cell(dataOut.totalNumStimuli,1);
    BinThresh = [];
    RespMatrixTraining = cell(length(selected_cells),1);
    RespMatrixTesting = cell(length(selected_cells),1);
    RespMatrixAllTrials = cell(length(selected_cells),1);
    
    %start calculations for each cell
    for c = 1:length(selected_cells)
        cell_N = selected_cells(c);
        
        %get only training trials 
        cell_N_resps = squeeze(response_matrix_all(:,cell_N,:));
        cell_N_resps_stimArray = cell(1,dataOut.totalNumStimuli);
        cell_N_resps_all = [];
        for stim = 1:dataOut.totalNumStimuli %for each stim
            cell_N_resps_stimArray{1,stim}=cell_N_resps(find(removedBlocks(:,stim)==0),stim);%put into a cell array all the responses that aren't a "badBlock" for each stim
            stim_trials_left = length(cell_N_resps_stimArray{1,stim});%how many trials are left (will be integer)
            stim_trials_leaveout = floor(stim_trials_left/folds);%find number of trials to leave out (will be integer)
            logic_vec = zeros(1,length(cell_N_resps_stimArray{1,stim}));
            logic_vec((test_set - 1)*stim_trials_leaveout+1:test_set*stim_trials_leaveout)=1;
            stim_trials_training{1,stim} = cell_N_resps_stimArray{1,stim}(~logic_vec);
            stim_trials_testing{1,stim} = cell_N_resps_stimArray{1,stim}(logical(logic_vec));
            cell_N_resps_all = [cell_N_resps_all stim_trials_training{1,stim}'];
        end
        RespMatrixTraining{c,1} = stim_trials_training;%store info for each cell
        RespMatrixTesting{c,1} = stim_trials_testing;
        RespMatrixAllTrials{c,1} = cell_N_resps_stimArray;
        
        %get sorted list of all responses in training data       
        cell_N_resps_sorted = sort(cell_N_resps_all); %sort all responses from training data 
        total_range = cell_N_resps_sorted(length(cell_N_resps_sorted)) - cell_N_resps_sorted(1); %get total range of responses
        bin_width = round(total_range/bins,3); %get the width of each bin in terms of response mag
        
        %get bin thresholds for this cell
        cell_N_bin_thresh = zeros(1,bins-1);
        for b = 1:bins-1
            bin = round((length(cell_N_resps_sorted)/bins)*b);
            cell_N_bin_thresh(b) = cell_N_resps_sorted(bin);
        end
        BinThresh = [BinThresh ; cell_N_bin_thresh];
        
        %start looping through each stim
        for stim = 1:dataOut.totalNumStimuli
            cell_N_stim_resp = RespMatrixTraining{c,1}{1,stim};
            cell_N_stim_resp_sorted = sort(cell_N_stim_resp);
            
            %get probability for each bin
            for i = 1:bins
                if i == 1
                    numbers = find(cell_N_stim_resp_sorted<cell_N_bin_thresh(i));%for the first bin, just find everything below the first bin thresh
                elseif i == bins
                    numbers = find(cell_N_stim_resp_sorted>cell_N_bin_thresh(i-1));%for the last bin, just find everything after the last bin tresh
                else
                    numbers = find(cell_N_stim_resp_sorted>=cell_N_bin_thresh(i-1) & cell_N_stim_resp_sorted<cell_N_bin_thresh(i));%for the rest of bins, find between that bin thresh and the next
                end
                %cell_N_prob = length(numbers)/length(cell_N_stim_resp);%calc probability for that bin
                cell_N_prob = (length(numbers)+.5)/(length(cell_N_stim_resp)+.5);%calc probability for that bin using new fix%%%%%%%%%%%%%%%
                ProbDist{stim,1}{c,i} = cell_N_prob;
            end               
        end
    end

end

