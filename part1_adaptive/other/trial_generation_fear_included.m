%coe part 1 adaptive
% Finnish data
%i = excerpts, j = false excerpts, k = emotions
%assumes emotion ratings from 7:11 column
emo_idxs = [7,8,9,10,11];
dist_ratings = struct;
data = readtable('data/FIN_merged_rawdata.csv');
format = readtable('trials.csv');
mean_scores = readtable('data/mean_ratings_set2.xls');
ratings = mean_scores{:,[5,6,7,8,9]};
fear_ratings = mean_scores{:,"fear"};
emo_labels = mean_scores.Properties.VariableNames([5,6,7,8,9]);
for i = 1:110
    excerpt_data = data(find(data{:,2} == i),:);
    mean_excerpt_values(i,:) = nanmean(excerpt_data{:,emo_idxs});
    [sorted,sorted_labels] = sort(mean_excerpt_values(i,:),'descend');
    [target_emotion,target_label_idx] = max(mean_excerpt_values(i,:));
    dists = target_emotion - sorted(2:end);
    fear_dists = fear_ratings(i)-mean_scores{i,[5,7,8,9]};
    totaldists(i) = mean(dists);
    emo_category(i) = target_label_idx;
    
    for j = 1:4
        difficulty = sum(excerpt_data{:,emo_idxs(sorted_labels(1))}>...
            excerpt_data{:,emo_idxs(sorted_labels(1+j))})/height(excerpt_data);
        false_label_difficulty = sum(excerpt_data{:,emo_idxs(sorted_labels(1))}<...
            excerpt_data{:,emo_idxs(sorted_labels(1+j))})/height(excerpt_data);
        [h,p,ci,stats] = ttest(excerpt_data{:,emo_idxs(sorted_labels(1))},...
            excerpt_data{:,emo_idxs(sorted_labels(1+j))});
        dist_ratings([i-1]*4+j).Label = strjoin([emo_labels(target_label_idx),emo_labels(sorted_labels(1+j))]);
        dist_ratings([i-1]*4+j).Name = i;
        dist_ratings([i-1]*4+j).TargetEmo = sorted(1);
        dist_ratings([i-1]*4+j).ComparisonEmo = sorted(1+j);
        dist_ratings([i-1]*4+j).Difficulty = difficulty;
        dist_ratings([i-1]*4+j).DifficultyInverted = false_label_difficulty;
        dist_ratings([i-1]*4+j).Distance = dists(j);
        dist_ratings([i-1]*4+j).Ttest = stats.tstat;
        dist_ratings([i-1]*4+j).Label1 = target_label_idx
        dist_ratings([i-1]*4+j).Label2 = sorted_labels(j+1)
        dist_ratings([i-1]*4+j).Soundtrack = mean_scores{i,'Soundtrack'};
    end
end
dist_ratings=struct2table(dist_ratings);
%remove fear trials
for i =1:height(dist_ratings)
    fear_idxs(i) = contains(dist_ratings.Label{i},'fear');
end
dist_ratings=dist_ratings(~fear_idxs,:);

%remove duplicate excerpts
[~,unique_excerpts] = unique(mean_scores{:,'IndexInSet1'});
unique_excerpts = sort(unique_excerpts);
excerpts_to_remove = [17,67,72,75,82,86,95,101];
idx = [];
for i =1:height(dist_ratings)
    idx(i) = ~sum(dist_ratings.Name(i)==excerpts_to_remove);
end
dist_ratings = dist_ratings(find(idx),:);

dist_ratings = sortrows(dist_ratings,'Ttest','descend');
%remove trials with low target emo means
p = prctile(dist_ratings.TargetEmo,10);
idx = find(dist_ratings.TargetEmo>p);
dist_ratings = dist_ratings(idx,:);

%create levels for each emotion separately
for k = 1:5
    cell_trials{k} = dist_ratings(find(dist_ratings.Label1 == k),:);
    levelLength = height(cell_trials{k})/16;
    for i = 1:height(cell_trials{k})
        cell_trials{k}.Level{i} = ceil(i/levelLength);
    end
end
trials=[];
for i = 1:5
    trials = [trials; cell_trials{i}];
end
trials = sortrows(trials,'Ttest','descend');

trials.Distance = trials.Ttest;

%figure,plot(sort(dist_ratings{:,'Distance'})), xlabel('Trials'),ylabel('Distance')
figure
subplot(1,3,1)
for i = 1:5
    hold on
    plot(cell_trials{i}.Difficulty)
end
legend(emo_labels),xlabel('Trials'),ylabel('Difficulty'),title('Percentages of true answer')
hold off
subplot(1,3,2)
for i = 1:5
    hold on
    plot(sort(cell_trials{i}.Distance,'descend'))
end
legend(emo_labels),xlabel('Trials'),ylabel('Distance'),title('Distances of true and false answer')
hold off
subplot(1,3,3)
for i = 1:5
    hold on
    plot(sort(cell_trials{i}.Ttest,'descend'))
end
legend(emo_labels),xlabel('Trials'),ylabel('Ttest value'),title('Ttest of true and false answer')
hold off

%compare Difficulty metrics
figure
subplot(1,3,1)
[rho,pval] = corr(trials.Difficulty,trials.Distance);
scatter(trials.Difficulty,rescale(trials.Distance,0,1))
xlabel('Percentage (percentage of true answer being higher)')
ylabel('Distance (distance of mean emotion ratings)')
title(['Scatterplot of Percentage and Distance: r = ', num2str(round(rho,2))])
subplot(1,3,2)
[rho,pval] = corr(trials.Ttest,trials.Distance);
scatter(trials.Difficulty,trials.Distance)
xlabel('Tstatistic ')
ylabel('Distance (distance of mean emotion ratings)')
title(['Scatterplot of Ttest and Distance: r = ', num2str(round(rho,2))])
subplot(1,3,3)
[rho,pval] = corr(trials.Difficulty,trials.Ttest);
scatter(trials.Difficulty,trials.Ttest)
xlabel('Percentage (percentage of true answer being higher)')
ylabel('Tstatistic')
title(['Scatterplot between Percentage and Ttest: r = ', num2str(round(rho,2))])

trials = removevars(trials,["Ttest","Difficulty","DifficultyInverted"]);
trials{:,'Trials'} = [1:height(trials)]';
%writetable(trials,'trials.csv')

%%find excerpts for part 2
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
%writematrix(excerpt,'Part2excerpts.csv')




