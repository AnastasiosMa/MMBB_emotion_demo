cd 'Downloads/adaptive_emotion_test'
data = readtable('mean_ratings_set2.xls');
ratings = data{:,[5,7,8,9]};
emo_labels = data.Properties.VariableNames([5,7,8,9]);
dist_ratings = struct;
for i = 1:length(ratings)
    [sorted,k] = sort(ratings(i,:),'descend');
    [target_emotion,target_label_idx] = max(ratings(i,:));
    dists = target_emotion - sorted(2:end);
    totaldists(i) = mean(dists);
    emo_categoty(i) = target_label_idx;
    for j = 1:3
    dist_ratings([i-1]*3+j).Label = strjoin([emo_labels(target_label_idx),emo_labels(k(1+j))]); 
    dist_ratings([i-1]*3+j).Name = data{i,'Number'};
    dist_ratings([i-1]*3+j).TargetEmo = target_emotion;
    dist_ratings([i-1]*3+j).ComparisonEmo = sorted(1+j);
    dist_ratings([i-1]*3+j).Distance = dists(j);
    dist_ratings([i-1]*3+j).Label1 = target_label_idx
    dist_ratings([i-1]*3+j).Label2 = k(j+1)
    dist_ratings([i-1]*3+j).Soundtrack = data{i,'Soundtrack'};
    end
end
dist_ratings=struct2table(dist_ratings);
dist_ratings = sortrows(dist_ratings,'Distance','descend');
dist_ratings = dist_ratings(1:end-1,:); 
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
trials = sortrows(trials,9,'ascend');

%figure,plot(sort(dist_ratings{:,'Distance'})), xlabel('Trials'),ylabel('Distance')
figure
for i = 1:4
    hold on
    plot(cell_trials{i}.Distance)
end
legend(emo_labels),xlabel('Trials'),ylabel('Distance')
hold off
writetable(trials,'trials.csv')

%%find excerpts for part 2
figure,plot(sort(totaldists)), xlabel('Trials'),ylabel('Distance'), title('Part2 distances')
low = prctile(totaldists,25);
mid = prctile(totaldists,50);
high = prctile(totaldists,75);

for i=1:4
    idx = find(emo_categoty==i);
    [~,tmp]=min(abs(totaldists(idx)-low));
    excerpt(i,1) = idx(tmp);
    [~,tmp]=min(abs(totaldists(idx)-mid));
    excerpt(i,2) = idx(tmp);
    [~,tmp]=min(abs(totaldists(idx)-high));
    excerpt(i,3) = idx(tmp);
end
writematrix(excerpt,'Part2excerpts.csv')
