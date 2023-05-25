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
        idx(k,i) = find(difficulty_scores(k,i)==trials.Distance,1);
    end
end
for i=1:height(trials)
    occurences(i) = sum(sum(idx==i));
end
figure
scatter(1-rescale(trials.Distance),occurences,'filled')