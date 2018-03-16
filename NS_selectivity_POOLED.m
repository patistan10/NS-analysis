%% Selectivity measure

expID = 'POOLED';
significance = 0.01;
thresholds = 200;
home=pwd;
cd ..
NS = load(sprintf('dataOut_NatScenes_%s.mat',expID));
G = load(sprintf('dataOut_gratings_%s.mat',expID));
cd(home);

stims = NS.dataOut.totalNumStimuli;
cells_picked = [1:NS.dataOut.totalNumCells];

%get matrix for responses and pvals for all responsive cells
resp_avg = NS.dataOut.stats.global.response_ACTUAL_avg_vals(cells_picked,:);
resp_pval_fdr = NS.dataOut.stats.global.response_average_pval_fdr(cells_picked,:);

AUC_all = zeros(length(cells_picked),1);
NS_selec = zeros(length(cells_picked),1);
%go through each cell and plot curve
for c = 1:length(cells_picked)
    n=cells_picked(c); %actual cell
    n_resps = resp_avg(n,:); %responses of that cell
    n_pval = resp_pval_fdr(n,:); %pvals for responses of that cell
    n_resps_sig = n_resps(find(n_pval<significance)); %only significant responses
    max_thresh = max(n_resps_sig); %max response
    n_resps_sig_norm = n_resps_sig./max_thresh; %responses normalized to max
    
    pass_thresh = zeros(thresholds+1,2);
    for t = 1:thresholds+1
        thresh = 1 - round(1/thresholds,4)*(t-1);
        pass_thresh(t,1) = thresh;%normalize to max response
        pass_thresh(t,2) = length(find(n_resps_sig_norm>=thresh))/stims;%as proportion of total num stims
    end 
    
    n_AUC = trapz(sort(pass_thresh(:,1)),pass_thresh(:,2)); %AUC proportion of stim with spacing of 1/thresholds
    AUC_all(n) = n_AUC;
    NS_selec(n) = 1-n_AUC;
    if isempty(n_resps_sig_norm)
        NS_selec(n) = -1;
    end
    
    plot(pass_thresh(:,2),pass_thresh(:,1))
    xlim([0,1])
    xlabel('proportion of stim with sig resp (p<.01)')
    ylabel('normalized mean response')
    title(sprintf('Cell %i , NS selectivity = %.2f',n,NS_selec(n)));
    saveas(gca,sprintf('NS_selec_%istim_cell%i.fig',stims,n))
    saveas(gca,sprintf('NS_selec_%istim_cell%i.png',stims,n))
    
    
%     %%%%%%test set%%%%%%
%     pass_thresh_test = pass_thresh;
%     pass_thresh_test(:,2) = [0:0.005:1];    
%     plot(pass_thresh_test(:,2),pass_thresh_test(:,1))
%     xlim([0,1])
%     AUC_test = trapz(pass_thresh_test(:,2));    
end

save('POOLED_NS_selectivity_data.mat','AUC_all','NS_selec')

%% get BW of all grating cells and compare with selectivity
cd ..
G = load('dataOut_Gratings_POOLED.mat');
Fits = load('gaussFit_results_POOLED.mat');
cd(home)

grating_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average;
grating_cells_bw = Fits.fit_oriBW(grating_cells)';
grating_cells_NSI = NS_selec(grating_cells);
% scatter(grating_cells_NSI,grating_cells_bw)
% ylim([0 90])

grating_cells_sharp = grating_cells(grating_cells_bw>4 & grating_cells_bw<14);
grating_cells_sharp_bw = Fits.fit_oriBW(grating_cells_sharp)';
grating_cells_sharp_NSI = NS_selec(grating_cells_sharp);
grating_cells_broad = grating_cells(grating_cells_bw>20 & grating_cells_bw<90);
grating_cells_broad_bw = Fits.fit_oriBW(grating_cells_broad)';
grating_cells_broad_NSI = NS_selec(grating_cells_broad);

save('broad_sharp_NSI.mat','grating_cells_sharp','grating_cells_sharp_bw','grating_cells_sharp_NSI',...
    'grating_cells_broad','grating_cells_broad_bw','grating_cells_broad_NSI');
% NS_cells = NS.dataOut.stats.global.responsive_cells_p001_fdr_average;
% NS_cells_NSI = NS_selec(NS_cells);

%% get cumulative graph
y_sharp = [1/length(grating_cells_sharp):1/length(grating_cells_sharp):1]';
sharp_NSI_sorted = sort(grating_cells_sharp_NSI);
y_broad = [1/length(grating_cells_broad):1/length(grating_cells_broad):1]';
broad_NSI_sorted = sort(grating_cells_broad_NSI);
%y_NS = [1/length(NS_cells):1/length(NS_cells):1]';
%NS_NSI_sorted = sort(NS_cells_NSI);

plot(sharp_NSI_sorted,y_sharp,'LineWidth',2);
hold on
plot(broad_NSI_sorted,y_broad,'LineWidth',2);
xlabel('NS selectivity (-1 means not responsive)');
ylabel('Proportion of population')
legend('sharp','broad','location','best')
title(sprintf('NS selectivity for Broad(%i,BW>20) and Sharp(%i,4>BW<14)',length(y_broad),length(y_sharp)))
%plot(NS_NSI_sorted,y_NS)
saveas(gca,'cumulative_dist_sharp_broad_NSI_POOLED.fig')
saveas(gca,'cumulative_dist_sharp_broad_NSI_POOLED.png')

