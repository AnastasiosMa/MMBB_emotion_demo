%coe part 1 adaptive
% Finnish data
%i = excerpts, j = false excerpts, k = emotions
%assumes emotion ratings from 7:11 column
emo_idxs = [7,8,9,10,11];
dist_ratings = struct;
data = readtable('../data/input/FIN_merged_rawdata.csv');
format = readtable('../data/output/trials.csv');
mean_scores = readtable('../data/input/mean_ratings_set2.xls');
ratings = mean_scores{:,[5,6,7,8,9]};
emo_labels = mean_scores.Properties.VariableNames([5,6,7,8,9]);
for i = 1:110
    excerpt_data = data(find(data{:,2} == i),:);
    mean_excerpt_values(i,:) = nanmean(excerpt_data{:,emo_idxs});
    [sorted,sorted_labels] = sort(mean_excerpt_values(i,:),'descend');
    [target_emotion,target_label_idx] = max(mean_excerpt_values(i,:));
    dists = target_emotion - sorted(2:end);
    totaldists(i) = mean(dists);
    emo_category(i) = target_label_idx;
end
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
%writematrix(excerpt,'../data/output/Part2excerpts.csv')




