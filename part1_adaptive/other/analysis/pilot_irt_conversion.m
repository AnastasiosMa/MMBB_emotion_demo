%publish('pilot_reliability.m','format','pdf','showCode',false);
%cd Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/analysis/
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
data = readtable('pilot_data.xlsx');
lvl_achieved_col = 121;
difficulty_col = 61;
level_ended = cell(1,2); difficulty_sum = cell(1,2); lvl_achieved = cell(1,2);
minimum_trial_num = 12;
trials = readtable('trials.csv');
emo_labels = {'anger','happy','sad','tender'};
%% Test trials information
disp(['Number of available trials: ', num2str(height(trials))])
%plot difficulty of trials for each emotion
figure
hold on
for i =1:4
    plot(sort(trials{trials.Label1==i,'Distance'}),'LineWidth',2);
    N_trials(i) = sum(trials.Label1==i);
end
xlabel('Trials');ylabel('Difficulty');title('Trial difficulty per emotion');
legend(emo_labels);

difficulty_scores = data{:,difficulty_col+1:difficulty_col+60};
for k=1:size(difficulty_scores,1)
    for i=1:size(difficulty_scores,2)
        trial_num(k,i) = find(difficulty_scores(k,i)==trials.Distance,1);
    end
end
for i=1:height(trials)
    occurences(i) = sum(sum(trial_num==i));
end
figure
scatter(1-rescale(trials.Distance),occurences,'filled')
%% Correctness ratio / Difficulty scatterplot
responses = data{:,2:61};
trial_num_v= trial_num(:); responses_v = responses(:);
correct_responses = trial_num_v(find(responses_v));
[counts_trial,group_trial] = groupcounts(sort(trial_num_v));
[counts_correct,group_correct] = groupcounts(sort(correct_responses));
for i = 1:height(trials)
    idx_trial = find(group_trial==i);
    idx_correct = find(group_correct==i);
    if idx_trial
        if idx_correct
            trial_response_ratio(i) = counts_correct(idx_correct)/counts_trial(idx_trial);
        else
            trial_response_ratio(i) = 0;
        end
    else
        trial_response_ratio(i) = nan;
    end
end

figure
scatter(1-rescale(trials.Distance),trial_response_ratio,'filled');