%% Create binary response matrix from Spanish and Finnish data
emo_labels = {'Angry','Fearful','Happy','Sad','Tender'};
excerpts_to_remove = [17,18,67,72,75,82,86,95,101];
exclude_trials_threshold = 10; %prctile to exclude trials
emo_N = length(emo_labels);
%% spanish
emo_idxs_spa = [15,16,14,18,17];
data_spa = readtable('../data/input/SPA_merged_rawdata.csv');
track_idx_spa = 22;
participant_spa = 3;

for i =1:height(data_spa)
   spa(i,:) = data_spa{i,[emo_idxs_spa, track_idx_spa,participant_spa]};
end
%change participant naming
spa(:,end) = spa(:,end)+200;
%% finnish
emo_idxs_fi = [7,8,9,10,11];
data_fi = readtable('../data/input/FIN_merged_rawdata.csv');
track_idx_fi = 2;
participant_fi = 6;

for i =1:height(data_fi)
   fi(i,:) = data_fi{i,[emo_idxs_fi, track_idx_fi,participant_fi]};
end 
%% Get Mean ratings
data = [fi;spa];
emo_idxs = 1:5;
for i = 1:110
    excerpt_data = data(data(:,end-1) == i,:);
    %remove nans
    excerpt_data = excerpt_data(~isnan(excerpt_data(:,1)),:);
    mean_excerpt_values(i,:) = nanmean(excerpt_data(:,emo_idxs));
    [B,I] = sort(mean_excerpt_values(i,:),'descend');
    target_emotion(i) = B(1);target_label_idx(i) = I(1);
    second_emotion(i) = B(2);second_label_idx(i) = I(2);
    third_emotion(i) = B(3);third_label_idx(i) = I(3);
end
remove_duplicates = mean_excerpt_values; 
remove_duplicates(excerpts_to_remove,:) = [];
remove_idxs = target_label_idx;
remove_idxs(excerpts_to_remove) = [];

%remove excerpts with low mean ratings 
%p = prctile(diag(remove_duplicates(:,remove_idxs)),exclude_trials_threshold);
%idx = find(diag(mean_excerpt_values(:,target_label_idx))>p);
%idx(find(sum(idx==excerpts_to_remove,2))) = [];
idx = [1:110]';
idx(excerpts_to_remove) = [];
%% Create binary matrix
binary_responses = nan(length(unique(data(:,end))),(emo_N-1)*length(idx));
%trial_name,participant_name
% Get percentage scores
subject_id = unique(data(:,end))';
for k = 1:length(subject_id)
    j=1;
    participant_data = data(find(data(:,end)==subject_id(k)),:);
    for i = idx'
        trial_data = participant_data(find(participant_data(:,end-1)==i),:);
        if ~isempty(trial_data) && ~any(isnan(trial_data(:,1)))
           wrong_emo_idx = setdiff(1:emo_N,target_label_idx(i));
           answers = trial_data(:,emo_idxs(target_label_idx(i)))>...
               trial_data(:,emo_idxs(wrong_emo_idx));
           binary_responses(k,(j-1)*(emo_N-1)+1:j*(emo_N-1)) = answers;
        end
        j = j+1;
    end
end
%% Create binary matrix of secondary emotion
subject_id = unique(data(:,end))';
secondary_binary_responses = nan(length(unique(data(:,end))),(emo_N-2)*length(idx));
for k = 1:length(subject_id)
    j=1;
    participant_data = data(find(data(:,end)==subject_id(k)),:);
    for i = idx'
        trial_data = participant_data(find(participant_data(:,end-1)==i),:);
        if ~isempty(trial_data) && ~any(isnan(trial_data(:,1)))
           wrong_emo_idx = setdiff(1:emo_N,[target_label_idx(i), second_label_idx(i)]);
           answers = trial_data(:,emo_idxs(second_label_idx(i)))>...
               trial_data(:,emo_idxs(wrong_emo_idx));
           secondary_binary_responses(k,(j-1)*(emo_N-2)+1:j*(emo_N-2)) = answers;
        end
        j = j+1;
    end
end
%% Create binary matrix of third emotion
subject_id = unique(data(:,end))';
third_binary_responses = nan(length(unique(data(:,end))),(emo_N-3)*length(idx));
for k = 1:length(subject_id)
    j=1;
    participant_data = data(find(data(:,end)==subject_id(k)),:);
    for i = idx'
        trial_data = participant_data(find(participant_data(:,end-1)==i),:);
        if ~isempty(trial_data) && ~any(isnan(trial_data(:,1)))
           wrong_emo_idx = setdiff(1:emo_N,[target_label_idx(i), second_label_idx(i),...
               third_label_idx(i)]);
           answers = trial_data(:,emo_idxs(third_label_idx(i)))>...
               trial_data(:,emo_idxs(wrong_emo_idx));
           third_binary_responses(k,(j-1)*(emo_N-3)+1:j*(emo_N-3)) = answers;
        end
        j = j+1;
    end
end
%% Plots
figure
plot(sort(nanmean(binary_responses)),'LineWidth',5)
ylabel('Response accuracy','FontSize',32);
xlabel('Items','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(nanmean(binary_responses))])
box on
grid on
%% Prepare data for exporting (primary mean items)
binary_responses = array2table(binary_responses);
%create trial info table
trial_name = 1:size(binary_responses,2);
%get target emotion
target_emo = repmat(target_label_idx(idx)',1,emo_N-1)';
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
%construct nationality vector
[participants{1,1:length(unique(fi(:,end)))}] = deal('fi');
[participants{1,end+1:end+length(unique(spa(:,end)))}] = deal('spa');

%construct gender and age vector
[~,ia] = unique(data_fi.id);
gender_fi = data_fi.Gender(ia);
age_fi = data_fi.Age(ia);
[~,ia] = unique(data_spa.Subject);
gender_spa = data_spa.Sex(ia);
gender_spa(strcmpi(gender_spa,'Mujer')) = {'Female'};
gender_spa(strcmpi(gender_spa,'Hombre')) = {'Male'};
gender = [gender_fi;gender_spa];

age = [age_fi; nan(size(binary_responses,1)-length(age_fi),1)+1];

%Add demographics to binary responses
binary_responses.Participant = participants';
binary_responses.Age = age;
binary_responses.Gender = gender;

table_idx = repmat(idx,1,emo_N-1)';
trial_info = table(labels(:),target_emo(:),incorrect_emo(:),...
    trial_name',table_idx(:),'VariableNames',{'Labels','TargetEmo',...
    'ComparisonEmo','TrialNum','Track110'});

%remove trials 0.5<x<1
d = nanmean(binary_responses{:,1:end-3});
i = find(d>0.5 & d<1);
binary_responses = [binary_responses(:,i), binary_responses(:,end-2:end)];
trial_info = trial_info(i,:);
%writetable(binary_responses,'../data/output/binary_responses/fear_binary_responses.csv');
%writetable(trial_info,'../data/output/binary_responses/fear_trial_info.csv');
%% Prepare data for exporting (SECONDARY mean items)
%Item Thresholds
upper_thx = 0.8;
min_mean_rating = 3;

secondary_binary_responses_t = array2table(secondary_binary_responses);
%create trial info table
second_trial_name = 1:size(secondary_binary_responses_t,2);
%get target emotion
second_emo = repmat(second_label_idx(idx)',1,emo_N-2)';
%get incorrect emotion
emo_options = [1:emo_N];
for i = 1:length(idx)
    second_incorrect_emo(:,i) = emo_options(emo_options~=target_emo(1,i) & ...
        emo_options~=second_emo(1,i));
end
for i = 1:length(idx)
    for k = 1:emo_N-2
        second_labels{k,i} = [emo_labels{second_emo(k,i)}, '-',...
            emo_labels{second_incorrect_emo(k,i)}];
    end
end

table_idx = repmat(idx,1,emo_N-2)';
second_trial_info = table(second_labels(:),second_emo(:),second_incorrect_emo(:),...
    second_trial_name',table_idx(:),'VariableNames',{'Labels','TargetEmo',...
    'ComparisonEmo','TrialNum','Track110'});

second_emotion_raiting = repmat(second_emotion(idx)',1,emo_N-2)';
%remove trials 0.75<x<1
d = nanmean(secondary_binary_responses_t{:,1:end});
i = find([d>0.5 & d<upper_thx] & [second_emotion_raiting(:)>min_mean_rating]');
secondary_binary_responses_t = secondary_binary_responses_t(:,i);
second_trial_info = second_trial_info(i,:);
%% Prepare data for exporting (THIRD mean items)
third_binary_responses_t = array2table(third_binary_responses);
%create trial info table
third_trial_name = 1:size(third_binary_responses_t,2);
%get target emotion
third_emo = repmat(third_label_idx(idx)',1,emo_N-3)';
%get incorrect emotion
emo_options = [1:emo_N];
for i = 1:length(idx)
    third_incorrect_emo(:,i) = emo_options(emo_options~=target_emo(1,i) & ...
        emo_options~=second_emo(1,i) & emo_options~=third_emo(1,i));
end
for i = 1:length(idx)
    for k = 1:emo_N-3
        third_labels{k,i} = [emo_labels{third_emo(k,i)}, '-',...
            emo_labels{third_incorrect_emo(k,i)}];
    end
end

table_idx = repmat(idx,1,emo_N-3)';
third_trial_info = table(third_labels(:),third_emo(:),third_incorrect_emo(:),...
    third_trial_name',table_idx(:),'VariableNames',{'Labels','TargetEmo',...
    'ComparisonEmo','TrialNum','Track110'});

third_emotion_raiting = repmat(third_emotion(idx)',1,emo_N-3)';
%remove trials 0.75<x<1
d = nanmean(third_binary_responses_t{:,1:end});
i = find([d>0.5 & d<upper_thx] & [third_emotion_raiting(:)>min_mean_rating]');
third_binary_responses_t = third_binary_responses_t(:,i);
third_trial_info = third_trial_info(i,:);
%% Secondary trial Plots
figure
hold on
plot(sort(nanmean(table2array(secondary_binary_responses_t)))','LineWidth',5)
plot(sort(nanmean(table2array(third_binary_responses_t)))','LineWidth',5)
ylabel('Response accuracy','FontSize',32);
xlabel('Items','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(nanmean(table2array(secondary_binary_responses_t)))])
title('Other Trials','Fontsize',32)
legend({'Second Highest','Third Highest'})
box on
grid on

figure
plot([sort(target_emotion)',sort(second_emotion)',sort(third_emotion)'],'LineWidth',5)
ylabel('Mean Rating','FontSize',32);
xlabel('Excerpts','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(target_emotion)])
title('Target and Secondary Mean Ratings','Fontsize',32)
legend({'Target','Second','Third'},'Location','best')
box on
grid on
%% Combine first and secondary trials
c_trial_info = [trial_info;second_trial_info; third_trial_info];
c_trial_info{:,'TrialNum'} = [1:height(c_trial_info)]';
c_binary_responses = [binary_responses(:,1:end-3) secondary_binary_responses_t,...
    third_binary_responses_t, binary_responses(:,end-2:end)];
%writetable(c_binary_responses,'../data/output/binary_responses/binary_responses.csv');
%writetable(c_trial_info,'../data/output/binary_responses/trial_info.csv');
%% Histograms
figure
subplot(1,2,1)
histogram(nanmean(table2array(binary_responses(:,1:end-3))));
ylabel('Count','FontSize',32);
xlabel('Response Accuracy','FontSize',24);
title(['Items of Intended Emotions N = ' num2str(size(binary_responses,2))],'FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
m = nanmean(nanmean(table2array(binary_responses(:,1:end-3))));
xline(m,'-',{['Mean Accuracy = ' num2str(round(m,2))]},'LineWidth',3);


subplot(1,2,2)
histogram(nanmean(table2array(c_binary_responses(:,1:end-3))));
ylabel('Count','FontSize',32);
xlabel('Response Accuracy','FontSize',24);
title(['All Items N = ' num2str(size(c_binary_responses,2))],'FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
m = nanmean(nanmean(table2array(c_binary_responses(:,1:end-3))));
xline(m,'-',{['Mean Accuracy = ' num2str(round(m,2))]},'LineWidth',3);