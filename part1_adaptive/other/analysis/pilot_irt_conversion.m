%publish('pilot_irt_conversion.m','format','pdf','showCode',false);
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
subplot(1,2,1)
scatter(1-rescale(trials.Distance),occurences,'filled')
xlabel('Normalised difficulty'); ylabel('Count')
title('Number of instances per trial')

%Occurences per level
for i =1:16
    occurences_per_level(i) = mean(occurences(find(trials.Level==i)));
end
subplot(1,2,2)
bar(occurences_per_level)
title('Mean number of occurence per trial on each level')
xlabel('Level'); ylabel('Mean count');
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
%fit linear model
[b,~,~,~,stats] = regress(trial_response_ratio',[1-rescale(trials.Distance),ones(height(trials),1)]);
p = polyfit(1-rescale(trials.Distance(~isnan(trial_response_ratio))),trial_response_ratio(~isnan(trial_response_ratio))',2);
figure
hold on
scatter(1-rescale(trials.Distance),trial_response_ratio,40,'filled');
xlabel('Difficulty');ylabel('Correct response ratio');
title({'Linear regression model',['R^2 = ' num2str(round(stats(1),2))]});
plot([0:1],b(1)*[0:1]+b(2),'LineWidth',2)
axis([0.2 1 0 1.2])
hold off