%% Compare Finnish and Spanish samples
mean_scores = readtable('../data/input/mean_ratings_set2.xls');
emo_labels = mean_scores.Properties.VariableNames([5,6,7,8,9]);
excerpts_to_remove = [17,18,67,72,75,82,86,95,101];
%% spanish
emo_idxs = [15,16,14,18,17];
data = readtable('../data/input/SPA_merged_rawdata.csv');
track_idx = 22;
for i = 1:110
    excerpt_data = data(find(data{:,track_idx} == i),:);
    %remove nans
    excerpt_data = excerpt_data(~isnan(excerpt_data{:,emo_idxs(1)}),:);
    mean_excerpt_values_spa(i,:) = nanmean(excerpt_data{:,emo_idxs});
end
mean_excerpt_values_spa(excerpts_to_remove,:) = [];
%% finnish
emo_idxs = [7,8,9,10,11];
data = readtable('../data/input/FIN_merged_rawdata.csv');
track_idx = 2;
for i = 1:110
    excerpt_data = data(find(data{:,track_idx} == i),:);
    mean_excerpt_values_fi(i,:) = nanmean(excerpt_data{:,emo_idxs});
end
mean_excerpt_values_fi(excerpts_to_remove,:) = [];
%% Correlate scores
rho = diag(corr(mean_excerpt_values_fi,mean_excerpt_values_spa));
figure
for i =1:4
    subplot(2,2,i)
    hold on
    scatter(mean_excerpt_values_fi(:,i),mean_excerpt_values_spa(:,i),'filled');
    plot([0,10],[0,10],'LineWidth',2)
    title({emo_labels{i},['Rho = ',num2str(round(rho(i),2))]});
    xlabel('Finnish mean scores');ylabel('Spanish mean scores');
    hold off
end