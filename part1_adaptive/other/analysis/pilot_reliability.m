%publish('pilot_reliability.m','format','pdf','showCode',false);
%cd Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/analysis/
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
data = readtable('pilot_data.xlsx');
lvl_achieved_col = 121;
difficulty_col = 61;
level_ended = cell(1,2); difficulty_sum = cell(1,2); lvl_achieved = cell(1,2);
minimum_trial_num = 15;
disp('3 Metrics used: Level ended, aggregate difficulties score, logistic regression slope')
%% Percentage of correct responses per level
responses = data{:,2:61};
current_level = data{:,122:181};
for trial = 1:60
    for k=1:size(data,1)
        if trial==1 || trial==31 && current_level(k,trial)==0
           current_level(k,trial)=1;
        elseif current_level(k,trial)==0
           current_level(k,trial) = current_level(k,trial-1);
        end
    end
end
disp(['Percentage of correct responses across all trials: ' num2str(round((sum(sum(responses))/numel(responses))*100,2))])
disp(['Percentage for 1st test: ' num2str(round((sum(sum(responses(:,1:30)))/numel(responses(:,1:30)))*100,2))])
disp(['Percentage for retest: ' num2str(round((sum(sum(responses(:,31:60)))/numel(responses(:,31:60)))*100,2))])

responses_v = responses(:);
current_level_v = current_level(:);
for i=1:16
    n_level_all(i) = sum(sum(current_level==i));
    correctness_all(i) = sum(responses_v(current_level_v==i))/n_level_all(i);
    n_level_test(i) = sum(sum(current_level(:,1:30)==i));
    correctness_test(i) = sum(responses_v(current_level_v==i))/n_level_test(i);
    n_level_retest(i) = sum(sum(current_level(:,31:60)==i));
    correctness_retest(i) = sum(responses_v(current_level_v==i))/n_level_retest(i);
end
figure
subplot(1,2,1)
plot([n_level_all;n_level_test;n_level_retest]','LineWidth',1)
xlabel('Level');ylabel('Count');title('Count of level instances')
subplot(1,2,2)
plot([correctness_all;correctness_test;correctness_retest]','LineWidth',1)
xlabel('Level');ylabel('Correct responses');title('Correct responses per level')
legend({'All trials','Test trials','Retest trials'})
%% Level ended
disp('The level the participant reached in the last trial')
%correct the 0's in level achieved data
for i =1:2
    lvl_achieved{i} = zeros(size(data,1),length(minimum_trial_num:30));
    for trial=1:30
        for k=1:size(data,1)
            if trial==1 && data{k,(i-1)*30+trial+lvl_achieved_col}==0
                lvl_achieved{i}(k,trial) = 1;
            elseif data{k,(i-1)*30+trial+lvl_achieved_col}==0
                jump_level= ceil(lvl_achieved{i}(k,trial-1)/(1/0.2));
                lvl_achieved{i}(k,trial) = lvl_achieved{i}(k,trial-1) - jump_level;
            else
                lvl_achieved{i}(k,trial) = data{k,(i-1)*30+trial+lvl_achieved_col};
            end
        end
    end
end
%get level the participant reached in the last trial
for i = 1:2
    level_ended{i} = zeros(size(data,1),length(minimum_trial_num:30));
    for trial=1:length(minimum_trial_num:30)
        for k=1:size(data,1)
            level_ended{i}(k,trial) = lvl_achieved{i}(k,trial+minimum_trial_num-1);
        end
    end
end
%Plot each participant's data
figure
for i = 1:2
    subplot(1,2,i)
    hold on
    plot(level_ended{i}','LineWidth',1)
    xlabel('Number of trials');ylabel('Level ended'); title('Level ended per participant')
    set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
    xtickangle(0)
    hold off
end
figure
for i = 1:2
    subplot(1,2,i)
    plot(median(level_ended{i}),'LineWidth',1)
    xlabel('Number of trials');ylabel('Level ended'); title('Median level ended per trial number')
    set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
    xtickangle(0)
    hold off
end
%% Aggregate difficulties Score
disp('Score based on difficulties of each trial. Score = sum(correct responses difficulty) - sum(1-incorrect responses difficulty)')
trial_data = readtable('trials.csv');
difficulty_scores = rescale(data{:,difficulty_col:difficulty_col+60});
difficulty_scores = 1-difficulty_scores;
for i = 1:2
    difficulty_sum{i} = zeros(size(data,1),length(minimum_trial_num:30));
    for trial=1:length(minimum_trial_num:30)
        for k=1:size(data,1)
            temp_correct=data{k,(i-1)*30+2:(i-1)*30+trial+minimum_trial_num-1};
            temp_difficulty=difficulty_scores(k,(i-1)*30+1:(i-1)*30+trial+minimum_trial_num-1);
            if any(temp_correct)
                difficulty_sum{i}(k,trial) = sum(temp_difficulty(find(temp_correct)));
            end
            if any(~temp_correct)
                difficulty_sum{i}(k,trial) = difficulty_sum{i}(k,trial)...
                    - sum(1-temp_difficulty(find(~temp_correct)));
            end
        end
    end
end
figure
for i = 1:2
    subplot(1,2,i)
    hold on
    plot(difficulty_sum{i}','LineWidth',1)
    xlabel('Number of trials');ylabel('Score'); title('Aggregate difficulty score per participant')
    set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
    xtickangle(0)
    hold off
end
%% Logistic regression slope
disp('Slope of logistic regression, X=trial difficulty, Y=correctness of response')
trial_data = readtable('trials.csv');
difficulty_scores = rescale(data{:,difficulty_col:difficulty_col+60});
difficulty_scores = 1-difficulty_scores;
for i = 1:2
    sigmoid_beta{i} = zeros(size(data,1),length(minimum_trial_num:30));
    for trial=1:length(minimum_trial_num:30)
        for k=1:size(data,1)
            temp_correct=data{k,(i-1)*30+2:(i-1)*30+trial+minimum_trial_num};
            temp_difficulty=difficulty_scores(k,(i-1)*30+1:(i-1)*30+trial+minimum_trial_num-1);
            b=mnrfit(temp_difficulty,temp_correct+1);
            sigmoid_beta{i}(k,trial) = b(2);
        end
    end
end
figure
for i = 1:2
    subplot(1,2,i)
    hold on
    plot(sigmoid_beta{i}')
    xlabel('Number of trials');ylabel('Score'); title('Logistic regression slope per participant')
    set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
    xtickangle(0)
    hold off
end

%% Test-retest reliability
labels = {'Level ended','Difficulty score','Logistic slope'};
disp('Calculating Cronbach alpha and test-retest correlation')
rho(:,1) = diag(corr(level_ended{1},level_ended{2}));
rho(:,2) = diag(corr(difficulty_sum{1},difficulty_sum{2}));
rho(:,3) = diag(corr(sigmoid_beta{1},sigmoid_beta{2}));

for i = 1:size(difficulty_sum{1},2)
    alpha(i,1) = cronbach([level_ended{1}(:,i),level_ended{2}(:,i)]);
    alpha(i,2) = cronbach([difficulty_sum{1}(:,i),difficulty_sum{2}(:,i)]);
    alpha(i,3) = cronbach([sigmoid_beta{1}(:,i),sigmoid_beta{2}(:,i)]);
end
figure
subplot(1,2,1)
plot(alpha,'LineWidth',1);
xlabel('Number of trials');ylabel('Alpha'); title('Cronbachs alpha per trial number')
set(gca,'XTick',1:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:30])
xtickangle(0)

subplot(1,2,2)
plot(rho,'LineWidth',1);
xlabel('Number of trials');ylabel('Correlation coefficient'); title('Correlation between tests per trial number')
set(gca,'XTick',1:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:30])
xtickangle(0)
legend(labels)

disp('Highest reliability score per method')
[max_alpha,idx_alpha] = max(alpha);
[max_rho,idx_rho] = max(rho);
trial_n = [minimum_trial_num:30];
res_alpha=array2table([max_alpha',trial_n(idx_alpha)'],'VariableNames',{'Cronbachs Alpha','Number of trials'},'RowNames',labels);
res_rho=array2table([max_rho',trial_n(idx_rho)'],'VariableNames',{'Correlation Coefficient','Number of trials'},'RowNames',labels);
disp(res_alpha);disp(res_rho);
