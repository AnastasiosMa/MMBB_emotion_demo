%cd Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/analysis/
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
pilot_data = readtable('pilot_data(with trials).csv');
new_trials = readtable('../data/output/binary_responses/final_trial_info.csv');
new_trials_difficulty = readtable('../data/output/binary_responses/irt_models/final_rasch_mirt.csv');
new_trials = [new_trials, new_trials_difficulty];
fname = 'pilot_trials.json';
fid = fopen(fname);
raw = fread(fid);
str = char(raw');
fclose(fid);
js= jsondecode(str);
old_trials = [];
for k  = 1:16
    for e = 1:4
        temp = struct2table(js.(['x',num2str(k)]).(['x',num2str(e)]));
        old_trials = [old_trials; temp];
    end
end
%% match old to new trials
emo_match = [1,3,4,5]; %emotion index of old emotions in the new emotions
for i = 1:height(old_trials)
    track_idx = find(new_trials.Track110==str2num(old_trials.Name{i}));
    target_emo_idx = track_idx(find(new_trials.TargetEmo(track_idx)==emo_match(str2num(old_trials.Label1{i}))));
    second_emo_idx = track_idx(find(new_trials.ComparisonEmo(track_idx)==emo_match(str2num(old_trials.Label2{i}))));
    if ~isempty(intersect(target_emo_idx,second_emo_idx))
        old2new(i,:) = [str2num(old_trials.Trials{i}),new_trials.TrialNum(intersect(target_emo_idx,second_emo_idx))];
    else
        old2new(i,:) = [str2num(old_trials.Trials{i}),nan];
    end
end
%% Recreate the binary responses table
emo_labels = {'Angry','Fearful','Happy','Sad','Tender'};
emo_N = length(emo_labels);
%spanish data
emo_idxs_spa = [15,16,14,18,17];
data_spa = readtable('../data/input/SPA_merged_rawdata.csv');
track_idx_spa = 22;
participant_spa = 3;
for i =1:height(data_spa)
    spa(i,:) = data_spa{i,[emo_idxs_spa, track_idx_spa,participant_spa]};
end
%change participant naming
spa(:,end) = spa(:,end)+200;
% finnish
emo_idxs_fi = [7,8,9,10,11];
data_fi = readtable('../data/input/FIN_merged_rawdata.csv');
track_idx_fi = 2;
participant_fi = 6;
for i =1:height(data_fi)
    fi(i,:) = data_fi{i,[emo_idxs_fi, track_idx_fi,participant_fi]};
end
% Get Mean ratings
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
idx = [1:110]';
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
% Create binary matrix of secondary emotion
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
% Create binary matrix of third emotion
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
%construct primary table
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
table_idx = repmat(idx,1,emo_N-1)';
trial_info = table(labels(:),target_emo(:),incorrect_emo(:),...
    trial_name',table_idx(:),'VariableNames',{'Labels','TargetEmo',...
    'ComparisonEmo','TrialNum','Track110'});
%construct secondary table
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
%third table
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
%combine trials
c_trial_info = [trial_info;second_trial_info; third_trial_info];
c_trial_info{:,'TrialNum'} = [1:height(c_trial_info)]';
c_binary_responses = [binary_responses, secondary_binary_responses_t,...
    third_binary_responses_t];
%% Calculate difficulty of missing old trials
missing_trials_idx = old2new(find(isnan(old2new(:,2))));
missing_trials = old_trials(missing_trials_idx,:);
missing_excerpts = unique(missing_trials.Name);
for i = 1:length(missing_excerpts)
    missing_excerpts{i} = str2num(missing_excerpts{i});
end
missing_excerpts = cell2mat(missing_excerpts);
%find indexes of missing excerpts in binary responses
emo_match = [1,3,4,5]; %emotion index of old emotions in the new emotions
for i = 1:height(missing_trials)
    track_idx = find(c_trial_info.Track110==str2num(missing_trials.Name{i}));
    target_emo_idx = track_idx(find(c_trial_info.TargetEmo(track_idx)==emo_match(str2num(missing_trials.Label1{i}))));
    second_emo_idx = track_idx(find(c_trial_info.ComparisonEmo(track_idx)==emo_match(str2num(missing_trials.Label2{i}))));
    missing_idx(i) = intersect(target_emo_idx,second_emo_idx);
end
missing_accuracy = nanmean(c_binary_responses{:,missing_idx});
%% Map missing old trials to new
new_trials_binary_responses = readtable('../data/output/binary_responses/final_binary_responses.csv');
new_trials_item_accuracy = nanmean(table2array(new_trials_binary_responses));
model = fitlm(new_trials_item_accuracy,new_trials{:,7});
y_pred = predict(model,missing_accuracy');
range_missing = max(y_pred)-min(y_pred);
range_new = max(new_trials{:,7})-min(new_trials{:,7});
range_ratio = range_new/range_missing;
y_pred = y_pred*range_ratio;

responses = table2array(pilot_data(:,2:61));
trials_taken = table2array(pilot_data(:,183:242));

%Add the missing trials in old2new
counter = 1;
for i = 1:size(old2new,1)
    if isnan(old2new(i,2))
        old2new(i,2) = length(new_trials{:,7})+counter;
        counter = counter+1;
    end
end

for k = 1:size(trials_taken,1)
    for i = 1:size(trials_taken,2)
        new_trial_number = old2new(find(trials_taken(k,i)==old2new(:,1)),2);
        if ~isempty(find(new_trials.TrialNum==new_trial_number))
            converted_trials_idx(k,i) = find(new_trials.TrialNum==new_trial_number);
        else
            converted_trials_idx(k,i) = new_trial_number;
        end
    end
end
item_difficulty_all = [new_trials{:,7};y_pred];
%% Create probability of correct and wrong sample answers
theta_step = 0.02;
theta_low = -6;
theta_high = 6;
theta_range = round(theta_low:theta_step:theta_high,2);
guessing = 0.5;
trialN = length(item_difficulty_all);
k=1;
for the = theta_range
    for j = 1:trialN
        p_correct(j,k) = guessing+(1-guessing)*(1/(1+exp(-(the-item_difficulty_all(j)))));
        p_star(j,k) = 1/(1+exp(-(the-item_difficulty_all(j))));
        p_incorrect(j,k) = 1-p_correct(j,k);
        product = p_correct(j,k)*p_incorrect(j,k);
        if k>1
            j_der1 = [p_correct(j,k)-p_correct(j,k-1)]./theta_step;
            information_test(j,k) = j_der1^2/product;
        end
    end
    k = k+1;
end
information_test(:,1) =  information_test(:,2);
%% calculate participant ability
for k = 1:size(responses,1)
    p_responses = responses(k,1:20);
    p_trials = converted_trials_idx(k,1:20);
    p_responses = p_responses(find(~isnan(p_trials)));
    p_trials = p_trials(find(~isnan(p_trials)));
    [th(k,1), th_est_idx,iter_N] = ml_optimizer(0,2,p_responses,p_trials,p_star,p_incorrect);
    p_responses = responses(k,31:50);
    p_trials = converted_trials_idx(k,31:50);
    p_responses = p_responses(find(~isnan(p_trials)));
    p_trials = p_trials(find(~isnan(p_trials)));
    [th(k,2), th_est_idx,iter_N] = ml_optimizer(0,2,p_responses,p_trials,p_star,p_incorrect);
end
corr(th(:,1),th(:,2))
[H,P,CI,STATS] = ttest(th(:,1),th(:,2));
%% Calculate cronbach
minimum_trial_num = 10;
for i = 1:2
    ability{i} = zeros(size(responses,1),length(minimum_trial_num:30));
    for trial=1:length(minimum_trial_num:30)
        for k=1:size(responses,1)
            p_responses = responses(k,(i-1)*30+1:(i-1)*30+trial+minimum_trial_num-1);
            p_trials = converted_trials_idx(k,(i-1)*30+1:(i-1)*30+trial+minimum_trial_num-1);
            [ability{i}(k,trial),~,~] = ml_optimizer(0,2,p_responses,p_trials,p_star,p_incorrect);
        end
    end
end
for i = 1:size(ability{1},2)
    alpha(i) = cronbach([ability{1}(:,i),ability{2}(:,i)]);
end

%calculate percentiles
centile = cell(1,2);
for i=1:2
    for k=1:size(ability{1},1)
        for trial = 1:size(ability{1},2)
            nless = sum(ability{i}(:,trial) < ability{i}(k,trial));
            nequal = sum(ability{i}(:,trial) == ability{i}(k,trial));
            centile{i}(k,trial) = 100 * (nless + 0.5*nequal) / length(ability{i}(:,trial));
            ability_deviation(i,trial) = nanmean(ability{i}(:,end)-ability{i}(:,trial));
            ability_deviation_std(i,trial) = nanstd(ability{i}(:,end)-ability{i}(:,trial));
        end
    end
end
rho_p = diag(corr(centile{1},centile{2}));
for i = 1:size(ability{1},2)
    alpha_p(i) = cronbach([centile{1}(:,i),centile{2}(:,i)]);
end
colors.blue = [0.1035, 0.4146, 0.6928];
colors.orange = [255, 127, 14] / 255;
figure
%subplot(1,2,1)
hold on
set(gca,'FontSize',32,'LineWidth',2)
p1 = plot(alpha_p','LineWidth',5,'Color',colors.blue);
p2 = plot(alpha','LineWidth',5,'Color',colors.orange);
xlabel('Number of Trials');
ylabel('Cronbachs alpha');
title('Test-retest Reliability')
set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
xtickangle(0)
xlim([1,30-minimum_trial_num+1])
grid on
box on
xline(20-minimum_trial_num+1,'-','LineWidth',5);
legend([p1,p2],{'Percentiles','Raw Scores'},'Location','best')
%text(-0.2,1,'a','Units','normalized','FontSize',30,'FontWeight','bold')
hold off

subplot(1,2,2)
hold on
set(gca,'FontSize',32,'LineWidth',2)
p1 = plot(ability_deviation(1,:)','LineWidth',5,'Color',colors.blue);
p2 = plot(ability_deviation(2,:)','LineWidth',5,'Color',colors.orange);
xlabel('Number of Trials');
ylabel('Mean Deviation (in θ)');
title('Ability Deviation')
set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
xtickangle(0)
xlim([1,30-minimum_trial_num+1])
ylim([-2,2])
grid on
box on
xline(20-minimum_trial_num+1,'-','LineWidth',5);
legend([p1,p2],{'Test','Retest'},'Location','best')
text(-0.2,1,'b','Units','normalized','FontSize',30,'FontWeight','bold')
hold off