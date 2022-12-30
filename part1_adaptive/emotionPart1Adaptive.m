cd 'Downloads/adaptive_emotion_test'
data = readtable('mean_ratings_set2.xls');
ratings = data{:,[5,7,8,9]};
emo_labels = data.Properties.VariableNames([5,7,8,9]);
dist_ratings = struct;
for i = 1:length(ratings)
    [sorted,k] = sort(ratings(i,:),'descend');
    [target_emotion,target_label_idx] = max(ratings(i,:));
    dists = target_emotion - sorted(2:end);
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
dist_ratings = sortrows(dist_ratings,'Distance');
levelLength = floor(height(dist_ratings)/15);
for k =1:height(dist_ratings)
        dist_ratings.Level{k} = floor(k/levelLength)+1;
end
%figure,plot(sort(dist_ratings{:,'Distance'}))
dist_ratings = dist_ratings(2:end,:); 
writetable(dist_ratings,'distancesPart1')


%% load distances
trials = readtable('distancesPart1');
N = 16;
emo = []; for k = 1:N/4, emo = [emo randperm(4)]; end
for k = 1:N
    temp = dist_ratings(find([cell2mat(dist_ratings.Level) == k & dist_ratings.Label1==emo(k)]),:);
    trial(k,:) = dist_ratings(randi(height(temp)),:);
end

