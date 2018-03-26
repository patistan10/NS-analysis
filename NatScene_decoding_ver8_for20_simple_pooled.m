function [AllFold_AllBins,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(ExpID,type,selected_cells,bins,folds)
% Natural Scene Decoding, DOESN'T SAVE OUTPUT
%   INPUTS:
%   ExpID = what you want the main save string to be
%   type = string name of group being used to decode (i.e. allResp, Broad)
%   selected_cells = vector of cells to use in this decoding
%   bins = number of bins to divide response into (EVEN DIVISION THRESHOLDING)
%   folds = number of folds for cross validation
%   OUTPUTS:
%   AllFold_AllBins = storage for output of each fold
%   selected_cells = cells that were used in this decoding
%   bins_accuracy = accuracy for each fold for given bin

for fold = 1:folds
    %get probability distributions
    [RespMatrix,RespMatrixTraining,RespMatrixTesting,RespMatrixAllTrials,ProbDist,BinThresh,removedBlocks_new,selected_cells]=get_prob_dist_ver8_for20_pooled(ExpID,selected_cells,bins,folds,fold);

    AllFold(1,fold).RespMatrix = RespMatrix;
    AllFold(1,fold).RespMatrixTraining = RespMatrixTraining;
    AllFold(1,fold).RespMatrixTesting = RespMatrixTesting;
    AllFold(1,fold).RespMatrixAllTrials = RespMatrixAllTrials;
    AllFold(1,fold).ProbDist = ProbDist;
    AllFold(1,fold).BinThresh = BinThresh;
    AllFold(1,fold).removedBlocks_new = removedBlocks_new;

    %get total number of test trials
    trials=[];%first get total number of test trials
    for q = 1:length(RespMatrixTesting{1, 1})
        trial = length(RespMatrixTesting{1,1}{1,q});
        trials = [trials trial];
    end
    total_test_trials = sum(trials);%total number of test trials

    %get which bin each neuron is for each test trial
    test_trial_bins = cell(1,length(RespMatrixTesting{1, 1}));%now will want each neuron's bin for each test trial for each stim
    for stim = 1:length(RespMatrixTesting{1, 1})
        total_test_trials = length(RespMatrixTesting{1,1}{1,stim});%get total number of test trials per stim
        for t = 1:total_test_trials
            for c = 1:length(selected_cells)
                cell_resp = RespMatrixTesting{c,1}{1,stim}(t);%given cell's response for a given stim for a given trial c = 1:length(selected_cells)
                cell_bin_thresh = BinThresh(c,:); %given cell's bin threshold
                for b = 1:bins
                    if b == 1 % if it's the first bin, only need to check if it's less than that bin
                        if cell_resp<=cell_bin_thresh(b)
                            cell_bin = b;
                        end
                    elseif b == bins %if it's the last bin, only need to check if greater than that bin
                        if cell_resp>cell_bin_thresh(b-1)
                            cell_bin = b;
                        end
                    else %for all other bins, need to check if greater than bin and less than next bin
                        if cell_resp>cell_bin_thresh(b-1) && cell_resp<=cell_bin_thresh(b)
                            cell_bin = b;
                        end
                    end
                end
                test_trial_bins{1,stim}{1,t}(c,1) = cell_bin;%store which bin the neuron was in
            end
        end
    end
    AllFold(1,fold).testTrialBins = test_trial_bins;

    %get probability of each stim
    trials=[];%first get total number of test trials
    for q = 1:length(RespMatrixAllTrials{1, 1})
        trial = length(RespMatrixAllTrials{1,1}{1,q});
        trials = [trials trial];
    end
    total_num_trials = sum(trials);%total number of trials
    stim_probs = [];
    for s = 1:length(RespMatrixAllTrials{1, 1}) %now get prob of each stim
        stim_p = length(RespMatrixAllTrials{1,1}{1,s})/total_num_trials;
        stim_probs = [stim_probs stim_p];%matrix w/ prob of each stim
    end

    %need to get each part, then will sum to get theta
    %equation: sum from i=1 to N of log(prob(ri|stim))+log(prob(stim))
    trial_argmax = cell(1,length(RespMatrixTesting{1, 1}));
    for test_trial_set = 1:length(RespMatrixTesting{1, 1})%this is just to store the test trial results in their corresponding stimulus cell{}
        total_test_trials = length(RespMatrixTesting{1,1}{1,test_trial_set});%get total number of test trials per stim
        for t = 1:total_test_trials
            for s = 1:length(RespMatrixTesting{1, 1})%now go through each stim so can use specific probability
                stim_prob = stim_probs(s);
                needs_summing = [];
                for c = 1:length(selected_cells)%for each cell
                    cell_N_bin = test_trial_bins{1,test_trial_set}{1,t}(c); % get the bin that the cell is in for this test trial
                    cell_N_resp_prob = ProbDist{s,1}{c,cell_N_bin};%now get the probability for that cell for that stim for that bin
                    %f = log(cell_N_resp_prob) + log(stim_prob);
                    f = log(cell_N_resp_prob);
                    needs_summing = [needs_summing f]; %store all of fs so can sum them
                end
                trial_argmax{1,test_trial_set}{1,t}(s,1) = sum(needs_summing)+ log(stim_prob);
            end
        end
    end
    AllFold(1,fold).Trial_Argmax = trial_argmax;

end
AllFold_AllBins{1,bins}=AllFold;

%get max theta of each trial for all the test trials
DecodingMatrix_All = [];
for fold = 1:folds
    n=0;
    DecodingMatrix = [];
    for stim = 1:length(AllFold(fold).RespMatrixTesting{1, 1})%this is to get which test trial belongs to which stim
        total_test_trials = length(AllFold(fold).RespMatrixTesting{1,1}{1,stim});%get total number of test trials per stim
        for t = 1:total_test_trials
            n=n+1;
            DecodingMatrix(n,1) = t; %first column will be the trial
            DecodingMatrix(n,2) = stim; %second column will be the actual stim
            [val ind]=max(AllFold(fold).Trial_Argmax{1, stim}{1, t}); %get the max value and the stim that belongs to that value 
            DecodingMatrix(n,3) = ind; %third column is predicted stim
            DecodingMatrix(n,4) = val; %fourth column is max theta for that trial
        end
    end
    DecodingMatrix_All = [DecodingMatrix_All ; DecodingMatrix];
end
AllFold_AllBins{2,bins} = DecodingMatrix_All;
%get accuracy for each fold
accuracy_all = [];
for fold = 1:folds
    set_size = length(DecodingMatrix_All)/folds;
    DM_chosen = DecodingMatrix_All([(fold-1)*set_size + 1:fold*(set_size)],:);
    accurate = [];
    for i = 1:length(DecodingMatrix_All)/folds
        real = DM_chosen(i,2);
        guess = DM_chosen(i,3);

        if real==guess
            accurate(i) = 1;
        else
            accurate(i) = 0;
        end
    end
    accuracy = sum(accurate)/set_size;
    accuracy_all = [accuracy_all accuracy];
    bins_accuracy{1,bins} = accuracy_all; %save accuracy for each bin size (bin size corresponds to location in cell array)
end

%fprintf('Bin %d finished\n',bins);


%save(sprintf('NBdecoding_%istim_n%i_%s.mat',size(RespMatrix,3),length(selected_cells),type),'AllFold_AllBins','selected_cells','bins_accuracy');

end