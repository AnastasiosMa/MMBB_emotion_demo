%% Combine Finnish and Spanish samples
mean_scores = readtable('../data/input/mean_ratings_set2.xls');
emo_labels = mean_scores.Properties.VariableNames([5,7,8,9]);
excerpts_to_remove = [17,18,67,72,75,82,86,95,101];
fear_ratings = mean_scores{:,"fear"};
exclude_trials_threshold = 10;
%% spanish
emo_idxs_spa = [15,14,18,17];
data_spa = readtable('../data/input/SPA_merged_rawdata.csv');
track_idx_spa = 22;

for i =1:height(data_spa)
   data{1}(i,:) =  [data_spa{i,emo_idxs_spa},data_spa{i,track_idx_spa}]; 
end
%% finnish
emo_idxs_fi = [7,9,10,11];
data_fi = readtable('../data/input/FIN_merged_rawdata.csv');
track_idx_fi = 2;

for i =1:height(data_fi)
   data{2}(i,:) =  [data_fi{i,emo_idxs_fi},data_fi{i,track_idx_fi}]; 
end

data = [data{1};data{2}];
%% Get difficulty scores
emo_idxs = [1:4];
dist_ratings = struct;
for i = 1:110
    excerpt_data = data(data(:,end) == i,:);
    %remove nans
    excerpt_data = excerpt_data(~isnan(excerpt_data(:,1)),:);
    mean_excerpt_values(i,:) = nanmean(excerpt_data(:,emo_idxs));
    [sorted,sorted_labels] = sort(mean_excerpt_values(i,:),'descend');
    [target_emotion,target_label_idx] = max(mean_excerpt_values(i,:));
    dists = target_emotion - sorted(2:end);
    fear_dists = fear_ratings(i)-mean_scores{i,[5,7,8,9]};
    if any(fear_dists<0) %if fear not highest rated emotion
       totaldists(i) = mean([dists fear_dists target_emotion-fear_ratings(i)]);
       emo_category(i) = target_label_idx;
    else %if fear highest rated emotion
       totaldists(i) = mean(fear_dists);
       emo_category(i) = 5; 
    end
    
    for j = 1:3
        percentage_correct = sum(excerpt_data(:,emo_idxs(sorted_labels(1)))>...
            excerpt_data(:,emo_idxs(sorted_labels(1+j))))/size(excerpt_data,1);
        percentage_false = sum(excerpt_data(:,emo_idxs(sorted_labels(1)))<...
            excerpt_data(:,emo_idxs(sorted_labels(1+j))))/size(excerpt_data,1);
        [h,p,ci,stats] = ttest(excerpt_data(:,emo_idxs(sorted_labels(1))),...
            excerpt_data(:,emo_idxs(sorted_labels(1+j))));
        dist_ratings([i-1]*3+j).Label = strjoin([emo_labels(target_label_idx),emo_labels(sorted_labels(1+j))]);
        dist_ratings([i-1]*3+j).Name = i;
        dist_ratings([i-1]*3+j).TargetEmo = sorted(1);
        dist_ratings([i-1]*3+j).ComparisonEmo = sorted(1+j);
        dist_ratings([i-1]*3+j).PercentageCorrect = percentage_correct;
        dist_ratings([i-1]*3+j).PercentageFalse = percentage_false;
        dist_ratings([i-1]*3+j).Distance = dists(j);
        dist_ratings([i-1]*3+j).Ttest = stats.tstat;
        dist_ratings([i-1]*3+j).Label1 = target_label_idx
        dist_ratings([i-1]*3+j).Label2 = sorted_labels(j+1)
        dist_ratings([i-1]*3+j).Soundtrack = mean_scores{i,'Soundtrack'};
    end
end
dist_ratings=struct2table(dist_ratings);
dist_ratings = sortrows(dist_ratings,'Ttest','descend');
%% Preprocess-remove trials
%remove duplicate excerpts
[~,unique_excerpts] = unique(mean_scores{:,'IndexInSet1'});
unique_excerpts = sort(unique_excerpts);
excerpts_to_remove = [17,67,72,75,82,86,95,101];
idx = [];
for i =1:height(dist_ratings)
    idx(i) = ~sum(dist_ratings.Name(i)==excerpts_to_remove);
end
dist_ratings = dist_ratings(find(idx),:);
%remove trials with low target emo means
p = prctile(dist_ratings.TargetEmo,exclude_trials_threshold);
idx = find(dist_ratings.TargetEmo>p);
%dist_ratings = dist_ratings(idx,:);
%% Create levels for each emotion 
for k = 1:4
    cell_trials{k} = dist_ratings(find(dist_ratings.Label1 == k),:);
    levelLength = height(cell_trials{k})/16;
    for i = 1:height(cell_trials{k})
        cell_trials{k}.Level{i} = ceil(i/levelLength);
    end
end
trials=[];
for i = 1:4
    trials = [trials; cell_trials{i}];
end
trials = sortrows(trials,'Ttest','descend');
%% Plots
%figure,plot(sort(dist_ratings{:,'Distance'})), xlabel('Trials'),ylabel('Distance')
figure
subplot(1,3,1)
for i = 1:4
    hold on
    plot(cell_trials{i}.PercentageCorrect)
end
legend(emo_labels),xlabel('Trials'),ylabel('Difficulty'),title('Percentages of true answer')
hold off
subplot(1,3,2)
for i = 1:4
    hold on
    plot(sort(cell_trials{i}.Distance,'descend'))
end
legend(emo_labels),xlabel('Trials'),ylabel('Distance'),title('Distances of true and false answer')
hold off
subplot(1,3,3)
for i = 1:4
    hold on
    plot(sort(cell_trials{i}.Ttest,'descend'))
end
legend(emo_labels),xlabel('Trials'),ylabel('Ttest value'),title('Ttest of true and false answer')
hold off

%compare Difficulty metrics
figure
subplot(1,3,1)
[rho,pval] = corr(trials.PercentageCorrect,trials.Distance);
scatter(trials.PercentageCorrect,rescale(trials.Distance,0,1),80,'filled')
xlabel('Response Accuracy','FontSize',26)
ylabel('Distance (mean emotion ratings)','FontSize',26)
title(['Response accuracy and Distance: r = ', num2str(round(rho,2))],...
    'FontSize',20)
set(gca,'FontSize',20,'LineWidth',3)
box on
grid on

subplot(1,3,2)
[rho,pval] = corr(trials.Ttest,trials.Distance);
scatter(trials.Ttest,trials.Distance,80,'filled')
xlabel('Tstatistic ','FontSize',26)
ylabel('Distance (mean emotion ratings)','FontSize',26)
title(['Ttest and Distance: r = ', num2str(round(rho,2))],...
    'FontSize',20)
set(gca,'FontSize',20,'LineWidth',3)
box on
grid on

subplot(1,3,3)
[rho,pval] = corr(trials.PercentageCorrect,trials.Ttest);
scatter(trials.PercentageCorrect,trials.Ttest,80,'filled')
xlabel('Response Accuracy','FontSize',26)
ylabel('Tstatistic','FontSize',26)
title(['Response accuracy and Ttest: r = ', num2str(round(rho,2))],...
    'FontSize',20)
set(gca,'FontSize',20,'LineWidth',3)
box on
grid on
%% Save trials
trials.Distance = trials.Ttest;
%trials = removevars(trials,["Ttest","PercentageCorrect","PercentageFalse"]);
trials{:,'Trials'} = [1:height(trials)]';
%writetable(trials,'../data/output/trials.csv')

%% Find excerpts for part 2
figure,plot(sort(totaldists)), xlabel('Trials'),ylabel('Distance'), title('Part2 distances')
low = prctile(totaldists,25);
mid = prctile(totaldists,50);
high = prctile(totaldists,75);

for i=1:length(unique(emo_category))
    idx = find(emo_category==i);
    [~,tmp]=min(abs(totaldists(idx)-low));
    excerpt(i,1) = idx(tmp);
    idx(tmp)=[];
    [~,tmp]=min(abs(totaldists(idx)-mid));
    excerpt(i,2) = idx(tmp);
    idx(tmp)=[];
    [~,tmp]=min(abs(totaldists(idx)-high));
    excerpt(i,3) = idx(tmp);
    idx(tmp)=[];
end
%writematrix(excerpt,'data/output/Part2excerpts.csv')
