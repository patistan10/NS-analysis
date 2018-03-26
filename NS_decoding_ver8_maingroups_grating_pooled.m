%% Decode using specific groups of neurons

expID = 'POOLED';
%gratingsdir = sprintf('H:/ProcessedDataArchive/Pati/NatScene2/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings',expID,expID);
%nsdir = sprintf('H:/ProcessedDataArchive/Pati/NatScene2/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes',expID,expID);
folds = 10;
%load hartley and OSI data
home = pwd;
cd ../..
gratings = load('dataOut_Gratings_POOLED.mat');
NS = load('dataOut_NatScenes_POOLED.mat');
cd(home)

%% get each group
grating_cells = gratings.dataOut.stats.global.responsive_cells_p001_fdr_average';
grating_index = gratings.dataOut.stats.global.responsive_cells_p001_fdr_average_index';
ns_cells = NS.dataOut.stats.global.responsive_cells_p001_fdr_average';
ns_index = NS.dataOut.stats.global.responsive_cells_p001_fdr_average_index';
all_resp_cells = find(grating_index==1 | ns_index==1);
all_cells = [1:NS.dataOut.totalNumCells];

group_grating = grating_cells;
group_nongrating = find(grating_index==0);
group_nsresp = ns_cells;
group_allresp = all_resp_cells;
group_all = all_cells;

save('all_groups.mat','group_grating','group_nongrating','group_nsresp','group_allresp','group_all');

%% Decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
types = {'grating','non-grating','NSresp','allResp','all'};%,'diverse_newthresh','allHartley','allNSresp'};
groups{1} = group_grating;
groups{2} = group_nongrating;
groups{3} = group_nsresp;
groups{4} = group_allresp;
groups{5} = group_all;

for g = 1:5
    %do decoding with specific group
    type = types{g};
    selected_cells = groups{g};

    [AllFold_AllBins,selected_cells,bins_accuracy] = NatScene_decoding_ver8_for20_pooled(expID,type,selected_cells,folds);
    accuracies_2bin{g} = bins_accuracy{1,2};
    %get confusion matrix
    chosenBin = 2;
    real_v_guessed = AllFold_AllBins{2,chosenBin}(:,2:3);
    real_v_guessed_sorted = sort(real_v_guessed,2);
    total_stim = size(AllFold_AllBins{1,2}(1).RespMatrix,3);
    confusion_matrix = zeros(20);
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
    saveas(gcf,sprintf('confusionMatrix_%istim_n%i_bin%i_%s_hot.fig',total_stim,length(selected_cells),chosenBin,type));
    saveas(gcf,sprintf('confusionMatrix_%istim_n%i_bin%i_%s_hot.png',total_stim,length(selected_cells),chosenBin,type));
    save(sprintf('confusionMatrix_%istim_n%i_bin%i_%s',total_stim,length(selected_cells),chosenBin,type),'confusion_matrix','chosenBin');
    close all
end

%plot groups against eachother
figure('Position',[100 200 800 600])
%colors = {'b','c','r','g'};
colors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880],[0.6350, 0.0780, 0.1840]};
for g = 1:5
    %scatter(repmat(g,1,folds),accuracies_2bin{g},60,'MarkerEdgeColor',[colors{g}])
    hold on
    %scatter(g,mean(accuracies_2bin{g}),70,'k','LineWidth',2)
    scatter(g,mean(accuracies_2bin{g}),70,'LineWidth',2,'MarkerEdgeColor',[colors{g}],'MarkerFaceColor',[colors{g}])
end
xlim([0 6])
ylim([0 1])
xticks([1:1:5])
% labels{1}=sprintf('grating(%i)',length(groups{1}));
% labels{2}=sprintf('non-grating(%i)',length(groups{2}));
% labels{3}=sprintf('NSresp(%i)',length(groups{1}));
% labels{4}=sprintf('all(%i)',length(groups{1}));
xticklabels(types)
text(1-.1,0.04,sprintf('%i',length(groups{1})),'FontSize',14)
text(2-.1,0.04,sprintf('%i',length(groups{2})),'FontSize',14)
text(3-.1,0.04,sprintf('%i',length(groups{3})),'FontSize',14)
text(4-.1,0.04,sprintf('%i',length(groups{4})),'FontSize',14)
text(5-.1,0.04,sprintf('%i',length(groups{5})),'FontSize',14)
ylabel('Accuracy')
set(gca,'FontSize',14)
title(sprintf('Decoding Accuracy (%d neurons total)',gratings.dataOut.totalNumCells))
saveas(gcf,sprintf('%s_NBdecoding_%istim_2bins_allgroup.fig',expID,size(confusion_matrix,2)));
saveas(gcf,sprintf('%s_NBdecoding_%istim_2bins_allgroup.png',expID,size(confusion_matrix,2)));


