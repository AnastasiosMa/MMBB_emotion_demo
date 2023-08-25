%% Create binary response matrix from Spanish and Finnish data
mean_scores = readtable('data/input/mean_ratings_set2.xls');
emo_labels = mean_scores.Properties.VariableNames([5,7,8,9]);
excerpts_to_remove = [17,18,67,72,75,82,86,95,101];
fear_ratings = mean_scores{:,"fear"};
exclude_trials_threshold = 10;
emo_N = length(emo_labels);
%% spanish
emo_idxs_spa = [15,14,18,17];
data_spa = readtable('data/input/SPA_merged_rawdata.csv');
track_idx_spa = 22;

for i =1:height(data_spa)
   data{1}(i,:) =  [data_spa{i,emo_idxs_spa},data_spa{i,track_idx_spa}]; 
end
%% finnish
emo_idxs_fi = [7,9,10,11];
data_fi = readtable('data/input/FIN_merged_rawdata.csv');
track_idx_fi = 2;

for i =1:height(data_fi)
   data{2}(i,:) =  [data_fi{i,emo_idxs_fi},data_fi{i,track_idx_fi}]; 
end

data = [data{1};data{2}];%% finnish
%% Get Mean ratings
emo_idxs = [1:4];

for i = 1:110
    excerpt_data = data(data(:,end) == i,:);
    %remove nans
    excerpt_data = excerpt_data(~isnan(excerpt_data(:,1)),:);
    mean_excerpt_values(i,:) = nanmean(excerpt_data(:,emo_idxs));
    [target_emotion(i),target_label_idx(i)] = max(mean_excerpt_values(i,:));
end

remove_duplicates = mean_excerpt_values; 
remove_duplicates(excerpts_to_remove,:) = [];
remove_idxs = target_label_idx;
remove_idxs(excerpts_to_remove) = [];

%remove excerpts with low mean ratings 
p = prctile(diag(remove_duplicates(:,remove_idxs)),exclude_trials_threshold);
idx = find(diag(mean_excerpt_values(:,target_label_idx))>p);
idx(find(sum(idx==excerpts_to_remove,2))) = [];
%% Create binary matrix Spanish
binary_responses_spa = nan(length(unique(data_spa.Subject)),(emo_N-1)*length(idx));
%trial_name,participant_name
% Get percentage scores
subject_id = unique(data_spa.Subject)';
for k = 1:length(subject_id)
    j=1;
    participant_data = data_spa(find(data_spa.Subject==subject_id(k)),:);
    for i = idx'
        trial_data = participant_data(find(participant_data.Track110==i),:);
        if ~isempty(trial_data) && ~any(isnan(trial_data{:,emo_idxs_spa}))
           wrong_emo_idx = setdiff(1:emo_N,target_label_idx(i));
           answers = trial_data{:,emo_idxs_spa(target_label_idx(i))}>...
               trial_data{:,emo_idxs_spa(wrong_emo_idx)};
           binary_responses_spa(k,(j-1)*(emo_N-1)+1:j*(emo_N-1)) = answers;
        end
        j = j+1;
    end
end
%% Create binary matrix Finnish
binary_responses_fi = nan(length(unique(data_fi.id)),(emo_N-1)*length(idx));
%trial_name,participant_name
% Get percentage scores
subject_id = unique(data_fi.id)';
for k = 1:length(subject_id)
    j=1;
    participant_data = data_fi(find(data_fi.id==subject_id(k)),:);
    for i = idx'
        trial_data = participant_data(find(participant_data.Track==i),:);
        if ~isempty(trial_data) && ~any(isnan(trial_data{:,emo_idxs_fi}))
           wrong_emo_idx = setdiff(1:emo_N,target_label_idx(i));
           answers = trial_data{:,emo_idxs_fi(target_label_idx(i))}>...
               trial_data{:,emo_idxs_fi(wrong_emo_idx)};
           binary_responses_fi(k,(j-1)*(emo_N-1)+1:j*(emo_N-1)) = answers;
        end
        j = j+1;
    end
end
%combine trials
binary_responses = [binary_responses_fi;binary_responses_spa];
%% Plots
figure
plot(sort(nanmean(binary_responses)),'LineWidth',5)
ylabel('Response accuracy','FontSize',24);
xlabel('Items','FontSize',24);
set(gca,'FontSize',24,'LineWidth',2)
xlim([1 length(nanmean(binary_responses))])
box on
grid on
%% Export data
binary_responses = array2table(binary_responses);
[participants{1,1:size(binary_responses_fi,1)}] = deal('fi');
[participants{1,end+1:end+size(binary_responses_spa,1)}] = deal('spa');

binary_responses.Participant = participants';

%get gender information
[~,ia] = unique(data_fi.id);
gender_fi = data_fi.Gender(ia);
[~,ia] = unique(data_spa.Subject);
gender_spa = data_spa.Sex(ia);
gender_spa(strcmpi(gender_spa,'Mujer')) = {'Female'};
gender_spa(strcmpi(gender_spa,'Hombre')) = {'Male'};
gender = [gender_fi;gender_spa];

%create trial info table
trial_name = 1:size(binary_responses,2);
%get target emotion
target_emo = repmat(target_label_idx(idx)',1,3)';
%get incorrect emotion
emo_options = [1:emo_N];
for i = 1:length(idx)
    incorrect_emo(:,i) = emo_options(emo_options~=target_emo(1,i));
end
for i = 1:length(idx)
    for k = 1:emo_N-1
        labels{k,i} = [emo_labels{target_emo(k,i)}, '-',...
            emo_labels{incorrect_emo(k,i)}];
    end
end
idx = repmat(idx,1,3)';
trial_info = table(labels(:),target_emo(:),incorrect_emo(:),...
    trial_name',idx(:),'VariableNames',{'Labels','TargetEmo',...
    'ComparisonEmo','TrialNum','Track110'});

%writetable(binary_responses,'data/output/binary_responses/binary_responses.csv');
%writetable(trial_info,'data/output/binary_responses/trial_info.csv');