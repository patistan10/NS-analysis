%% Decode using specific groups of neurons
expID = 'POOLED';
folds = 10;
bins = 3;
%load NatScene data
home = pwd;
cd ../..
load('gaussFit_results_POOLED.mat');
load('dataOut_Gratings_POOLED.mat');
cd(home)

grating_cells = dataOut.stats.global.responsive_cells_p001_fdr_average';

%% get each group
grouping = zeros(1,length(gaussFitResults)); %1 is S, 2 is M, 3 is B
thresh1=14;
thresh2=20;
for n=1:length(gaussFitResults)
    if fit_oriBW(n)>=4 && fit_oriBW(n)<thresh1
        grouping(n)=1; %sharp
    elseif fit_oriBW(n)>thresh2 && fit_oriBW(n)<=90
        grouping(n)=3; %broad
    elseif fit_oriBW(n)>=thresh1 && fit_oriBW(n)<=thresh2
        grouping(n)=2; %medium
    end
end
grouping = grouping(grating_cells);%only want cells that have sig resp
grating_bw = fit_oriBW(grating_cells);%bw for cells with sig resp
%up to now just looking at grating cells, but want to also get rid of ones
%with weird BW
grating_bw(grouping==0) = [];%get rid of cells where bw was <4 or >90
grating_cells_new = grating_cells;%get rid of cells where bw was <4 or >90
grating_cells_new(grouping==0) = [];%get rid of cells where bw was <4 or >90
grouping(grouping==0) = []; %need this for next step

cells_sharp = grating_cells_new(grouping==1); %all sharp cells
cells_med = grating_cells_new(grouping==2); %all medium cells
cells_broad = grating_cells_new(grouping==3); %all broad cells

%make sharp and broad group same size and as large as possible
target_group_size = min([length(cells_broad) length(cells_med) length(cells_sharp)]); %get group size
max_group_size = max([length(cells_broad) length(cells_med) length(cells_sharp)]);    %also get size of larger group

%%%% do decoding for multiples of larger group
%get which groups are larger than group size, will do multiple repeats
large_group_ind = zeros(1,3);
if length(cells_broad) > target_group_size
    large_group_ind(1) = 1;
end
if length(cells_med)> target_group_size
    large_group_ind(2) = 1;
end
if length(cells_sharp)> target_group_size
    large_group_ind(3) = 1;
end
types = ({'broad','medium','sharp'});
%rand_vec_rep = repmat(rand_vec,1,2);

for g=1:3    
    if large_group_ind(g)==1
        type = types{g};
        if g==1 %for broad
            rand_vec = randsample(length(cells_broad),length(cells_broad))';
            rand_vec_rep = repmat(rand_vec,1,2);
            for i = 1:30            
                vec_picked = rand_vec_rep(i:i+target_group_size - 1);            
                selected_cells = cells_broad(vec_picked);
                broad_groups_all{i,1} = selected_cells;
                [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,selected_cells,bins,folds);
                mean_acc = mean(bins_accuracy{1,bins});
                broad_groups_all{i,2} = mean_acc;
                fprintf('rep %d finished\n',i);
            end
            accs = cell2mat(broad_groups_all(:,2));
            m = mean(accs);
            diffs = abs(accs - m);
            [min_diff,ind] = min(diffs);
            group_broad = broad_groups_all{ind,1};
        elseif g==2 %for med
            rand_vec = randsample(length(cells_med),length(cells_med))';
            rand_vec_rep = repmat(rand_vec,1,2);
            for i = 1:30            
                vec_picked = rand_vec_rep(i:i+target_group_size - 1);            
                selected_cells = cells_med(vec_picked);
                med_groups_all{i,1} = selected_cells;
                [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,selected_cells,bins,folds);
                mean_acc = mean(bins_accuracy{1,bins});
                med_groups_all{i,2} = mean_acc;
                fprintf('rep %d finished\n',i);
            end
            accs = cell2mat(med_groups_all(:,2));
            m = mean(accs);
            diffs = abs(accs - m);
            [min_diff,ind] = min(diffs);
            group_med = med_groups_all{ind,1};
        elseif g==3 %for sharp
            rand_vec = randsample(length(cells_sharp),length(cells_sharp))';
            rand_vec_rep = repmat(rand_vec,1,2);
            for i = 1:30            
                vec_picked = rand_vec_rep(i:i+target_group_size - 1);            
                selected_cells = cells_sharp(vec_picked);
                sharp_groups_all{i,1} = selected_cells;
                [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,selected_cells,bins,folds);
                mean_acc = mean(bins_accuracy{1,bins});
                sharp_groups_all{i,2} = mean_acc;
                fprintf('rep %d finished\n',i);
            end
            accs = cell2mat(sharp_groups_all(:,2));
            m = mean(accs);
            diffs = abs(accs - m);
            [min_diff,ind] = min(diffs);
            group_sharp = sharp_groups_all{ind,1};
        end
    else %if not more than one option
        type = types{g};
        if g==1 %for broad
            group_broad = cells_broad;
            broad_groups_all{1,1} = group_broad;
            [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,group_broad,bins,folds);
            broad_groups_all{1,2} = mean(bins_accuracy{1,bins});
        elseif g==2 %for med
            group_med = cells_med;
            med_groups_all{1,1} = group_med;
            [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,group_med,bins,folds);
            med_groups_all{1,2} = mean(bins_accuracy{1,bins});
        elseif g==3 %for sharp
            group_sharp = cells_sharp;
            sharp_groups_all{1,1} = group_sharp;
            [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,group_sharp,bins,folds);
            sharp_groups_all{1,2} = mean(bins_accuracy{1,bins});
        end
    end
end

save('group_avg_broad_med_sharp.mat','sharp_groups_all','med_groups_all','broad_groups_all','group_broad','group_med','group_sharp');

%% get diverse group 

%first put cells in group based on OSI
limits = [0:5:90];
BW_bins = [];
Ori_bw_all = grating_bw;
for lim = 1: length(limits)-1
    cells_BW = grating_cells_new(find(Ori_bw_all>limits(lim) & Ori_bw_all<limits(lim+1)));
    BW_bins(lim) = length(cells_BW);
    cells_bins{1,lim} = cells_BW;
end


%get how many cells to grab from each OSI group
BW_bins_session = round(BW_bins./length(Ori_bw_all).*target_group_size);

%make sure will have same size group as other two
while sum(BW_bins_session)~=target_group_size
    warning('diverse group size does not match min bin size, fixing to match');
    if sum(BW_bins_session)> target_group_size
        [bins_max,ind]=max(BW_bins);
        BW_bins_session(ind) = BW_bins_session(ind)-1;
    else
        [bins_max,ind]=max(BW_bins);
        BW_bins_session(ind) = BW_bins_session(ind)+1;
    end
end

%get random cells from each group to make up distribution
%group_diverse = [];
for r = 1:30
    selected_cells = [];
    for i = 1: length(BW_bins_session)
        take = BW_bins_session(i);
        rand=randsample(BW_bins(i),BW_bins(i));
        selected_cells = [selected_cells cells_bins{1,i}(rand(1:take))];
    end
    %now deal with making sure the groups aren't the same
    if r>1
        for check = 1:r-1
           if isequal(diverse_groups_all{check,1},selected_cells)%if the groups match completely
               warning('group repeated, trying again');
              %repeat the process of choosing a new group
              for i = 1: length(BW_bins_session)
                take = BW_bins_session(i);
                rand=randsample(BW_bins(i),BW_bins(i));
                selected_cells = [selected_cells cells_bins{1,i}(rand(1:take))];
              end
              %also make chek = 1 so it goes through and checks it all
              %again
              check = 1;           
           end
        end
    end
    [~,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_simple_pooled(expID,type,selected_cells,bins,folds);
    diverse_groups_all{r,1} = selected_cells;
    mean_acc = mean(bins_accuracy{1,bins});
    diverse_groups_all{r,2} = mean_acc;
    fprintf('rep %d finished\n',r);
end
accs = cell2mat(diverse_groups_all(:,2));
m = mean(accs);
diffs = abs(accs - m);
[min_diff,ind] = min(diffs);
group_diverse = diverse_groups_all{ind,1};        

save('diverse_avg.mat','diverse_groups_all','group_diverse');

%plot groups against eachother using averages
mean_b = mean(cell2mat(broad_groups_all(:,2)));
mean_m = mean(cell2mat(med_groups_all(:,2)));
mean_s = mean(cell2mat(sharp_groups_all(:,2)));
mean_d = mean(cell2mat(diverse_groups_all(:,2)));
figure('Position',[100 200 800 600])
hold on
if size(broad_groups_all,1)>1
    scatter(repmat([1],1,30),cell2mat(broad_groups_all(:,2)))
end
if size(med_groups_all,1)>1
    scatter(repmat([2],1,30),cell2mat(med_groups_all(:,2)))
end
if size(sharp_groups_all,1)>1
    scatter(repmat([3],1,30),cell2mat(sharp_groups_all(:,2)))
end
scatter(repmat([4],1,30),cell2mat(diverse_groups_all(:,2)))
scatter(1,mean_b,80,'b','filled','LineWidth',2)
scatter(2,mean_m,80,'c','filled','LineWidth',2)
scatter(3,mean_s,80,'r','filled','LineWidth',2)
scatter(4,mean_d,80,'g','filled','LineWidth',2)
ylim([0 1])
xlim([0 5])
xticks([1:1:4])
xticklabels({'broad','medium','sharp','diverse'})
ylabel('Accuracy')
set(gca,'FontSize',16)
title(sprintf('Decoding Accuracy (%d neurons per group)',length(group_sharp)))
saveas(gcf,sprintf('%s_NBdecoding_%istim_n%i_%ibins_allgroup_avg.fig',expID,20,length(group_sharp),bins));
saveas(gcf,sprintf('%s_NBdecoding_%istim_n%i_%ibins_allgroup_avg.png',expID,20,length(group_sharp),bins));

%% Decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('specific_groups','dir')
    mkdir('specific_groups')
end

types = {'broad','medium','sharp','diverse'};
groups{1} = group_broad;
groups{2} = group_med;
groups{3} = group_sharp;
groups{4} = group_diverse;

for g = 1:length(groups)
    %do decoding with specific group
    type = types{g};
    selected_cells = groups{g};

    [AllFold_AllBins,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_pooled(expID,type,selected_cells,bins,folds);
    accuracies_bin{g} = bins_accuracy{1,bins};
    
    %get confusion matrix    
    real_v_guessed = AllFold_AllBins{2,bins}(:,2:3);
    real_v_guessed_sorted = sort(real_v_guessed,2);
    total_stim = size(AllFold_AllBins{1,bins}(1).RespMatrix,3);
    confusion_matrix = zeros(total_stim);
    for i = 1:total_stim
        stimnum = real_v_guessed(find(real_v_guessed(:,1)== i),:);
        for n = 1:total_stim
            stimguessed = find(stimnum(:,2)==n);
            confused_perc = length(stimguessed)/length(stimnum);
            confusion_matrix(n,i) = confused_perc;
        end
    end
    figure('Position',[100 200 1000 800])
    imagesc(confusion_matrix)
    h=colorbar;
    ylabel(h, 'P(PS|S)')
    colormap hot
    caxis([0 1])
    xlabel('True Stim (S)')
    ylabel('Predicted Stim (PS)')
    set(gca,'FontSize',16)
    title(sprintf('Confusion Matrix(%d neurons, %s)',length(selected_cells),type))
    saveas(gcf,sprintf('confusionMatrix_%istim_n%i_bin%i_%s.fig',total_stim,length(selected_cells),bins,type));
    saveas(gcf,sprintf('confusionMatrix_%istim_n%i_bin%i_%s.png',total_stim,length(selected_cells),bins,type));
    save(sprintf('confusionMatrix_%istim_n%i_bin%i_%s',total_stim,length(selected_cells),bins,type),'confusion_matrix','bins');
    close all
end

%plot groups against eachother
figure('Position',[100 200 800 600])
colors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880],[0.6350, 0.0780, 0.1840]};
for g = 1:length(groups)
    scatter(repmat(g,1,folds),accuracies_bin{g},60,'MarkerEdgeColor',[colors{g}])
    hold on
    %scatter(g,mean(accuracies_2bin{g}),70,'k','LineWidth',2)
    scatter(g,mean(accuracies_bin{g}),70,'LineWidth',2,'MarkerEdgeColor',[colors{g}],'MarkerFaceColor',[colors{g}])
end
xlim([0 5])
xticks([1:1:4])
xticklabels({'broad','medium','sharp','diverse'})
ylabel('Accuracy')
set(gca,'FontSize',16)
title(sprintf('Decoding Accuracy (%d neurons per group)',length(selected_cells)))
saveas(gcf,sprintf('NBdecoding_%istim_n%i_%ibins_allgroup.fig',size(confusion_matrix,2),length(selected_cells),bins));
saveas(gcf,sprintf('NBdecoding_%istim_n%i_%ibins_allgroup.png',size(confusion_matrix,2),length(selected_cells),bins));


