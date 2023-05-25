%publish('pilot_reliability.m','format','pdf','showCode',false);
%cd Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/analysis/
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
data = readtable('pilot_data.xlsx');
lvl_achieved_col = 121;
difficulty_col = 61;
level_ended = cell(1,2); difficulty_sum = cell(1,2); lvl_achieved = cell(1,2);
minimum_trial_num = 12;

%% Analysis Goals
disp(['1)Find the best metric to assess participants performance. 2)Estimate test-retest reliability.'...
    ,' 3)Find the optimal number of test trials based on reliability scores.'])
disp('Metrics used to assess participants performance: level ended, aggregate difficulties score, logistic regression slope')
disp('Reliability metrics: Cronbachs alpha, correlation')
disp(['Range of trials tested: ', num2str(minimum_trial_num), ' to 30'])
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

responses_test = responses(:,1:30);
responses_retest = responses(:,31:60);
current_level_test = current_level(:,1:30);
current_level_retest = current_level(:,31:60);
for i=1:16
    n_level_all(i) = sum(sum(current_level==i));
    correctness_all(i) = sum(sum(responses(current_level==i)))/n_level_all(i);
    n_level_test(i) = sum(sum(current_level_test==i));
    correctness_test(i) = sum(sum(responses_test(current_level_test==i)))/n_level_test(i);
    n_level_retest(i) = sum(sum(current_level_retest==i));
    correctness_retest(i) = sum(sum(responses_retest(current_level_retest==i)))/n_level_retest(i);
end
%interpolate nans
X = ~isnan(correctness_all);
Y = cumsum(X-diff([1,X])/2);
correctness_all = interp1(1:nnz(X),correctness_all(X),Y);
X = ~isnan(correctness_test);
Y = cumsum(X-diff([1,X])/2);
correctness_test = interp1(1:nnz(X),correctness_test(X),Y);
X = ~isnan(correctness_retest);
Y = cumsum(X-diff([1,X])/2);
correctness_retest = interp1(1:nnz(X),correctness_retest(X),Y);

figure
subplot(1,2,1)
plot([n_level_test;n_level_retest]','LineWidth',1)
xlabel('Level');ylabel('Count');title('Count of level instances')
subplot(1,2,2)
plot([correctness_test;correctness_retest]','LineWidth',1)
xlabel('Level');ylabel('Correct responses');title('Correct responses per level')
legend({'Test trials','Retest trials'})
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
disp('Score based on difficulty of each trial. Difficulties rescaled to [0,1] interval. Score = sum(correct responses difficulty) - sum(1-incorrect responses difficulty)')
trial_data = readtable('trials.csv');
difficulty_scores = rescale(data{:,difficulty_col+1:difficulty_col+60});
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
            difficulty_sum{i}(k,trial) = difficulty_sum{i}(k,trial)/(trial+minimum_trial_num);
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
difficulty_scores = rescale(data{:,difficulty_col+1:difficulty_col+60});
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
    xtickangle(0); ylim([0 10]);
    hold off
end
%plot outliers in logistic slopes
figure
hold on
b=mnrfit(difficulty_scores(17,1:15),data{17,2:16}+1);
x=linspace(0,1,100);
plot(x,1./(1+exp(+b(1)+b(2)*x)),'linewidth',2)
scatter(difficulty_scores(17,1:15),data{17,2:16},'k','filled')
title('Data of outlier slope (Participant 17, b = 1206) to understand the extreme values');xlabel('Difficulty');ylabel('Response')
ylim([0,1.5])
legend('Sigmoid')
hold off
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
set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
xtickangle(0)

subplot(1,2,2)
plot(rho,'LineWidth',1);
xlabel('Number of trials');ylabel('Correlation coefficient'); title('Correlation between tests per trial number')
set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
xtickangle(0)
legend(labels)

disp('Highest reliability score per method')
[max_alpha,idx_alpha] = max(alpha);
[max_rho,idx_rho] = max(rho);
trial_n = [minimum_trial_num:30];
res_alpha=array2table([max_alpha',trial_n(idx_alpha)'],'VariableNames',{'Cronbachs Alpha','Number of trials'},'RowNames',labels);
res_rho=array2table([max_rho',trial_n(idx_rho)'],'VariableNames',{'Correlation Coefficient','Number of trials'},'RowNames',labels);
disp(res_alpha);disp(res_rho);

%% Explore test-retest differences
disp('Paired sample ttests between test and retest aggregate scores. Positive tstatistics denote higher test scores,negative values higher retest scores.')
for i=1:size(difficulty_sum{1},2)
    [h,p,ci,stats] = ttest(difficulty_sum{1}(:,i),difficulty_sum{2}(:,i));
    ttests(i) = stats.tstat;
end
figure
plot(ttests,'LineWidth',1);
title('Tstatistics across different trial numbers')
xlabel('Number of trials');ylabel('Tstatistic');
set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
%plot difficulty scores between most similar and dissimilar trial number
%in terms of tstatistic
[sortedt,idxs] = sort(abs(ttests),'descend');
idxs = [idxs(1), idxs(end)];
%slope=difficulty_sum{1}(:,8)'/difficulty_sum{2}(:,8)'
figure
for i =1:2
subplot(1,2,i)
hold on
scatter(difficulty_sum{1}(:,idxs(i))',difficulty_sum{2}(:,idxs(i))','filled')
axis([0 1 0 1])
lsline()
%plot([7:13],slope*[7:13])
plot([0,1],[0,1])
hold off
legend({'','Least squares line','Diagonal'})
title(['Scores at ',...
    num2str(idxs(i)+minimum_trial_num-1), ' trials, tstat = ',num2str(round(ttests(idxs(i)),2))])
xlabel('Test Score');ylabel('Retest Score');
end
sgtitle('Participant scores for trials with largest and smallest t-statistics')

%% Percentile reliability 
disp('Participants aggregate scores are converted to percentiles (calculated separately for test and retest).')
% Compute centile
centile = cell(1,2);
for i=1:2
    for k=1:size(difficulty_sum{1},1)
        for trial = 1:size(difficulty_sum{1},2)
            nless = sum(difficulty_sum{i}(:,trial) < difficulty_sum{i}(k,trial));
            nequal = sum(difficulty_sum{i}(:,trial) == difficulty_sum{i}(k,trial));
            centile{i}(k,trial) = 100 * (nless + 0.5*nequal) / length(difficulty_sum{i}(:,trial));
        end
    end
end
rho_p = diag(corr(centile{1},centile{2}));
for i = 1:size(difficulty_sum{1},2)
    alpha_p(i) = cronbach([centile{1}(:,i),centile{2}(:,i)]);
end
figure
plot([alpha_p',alpha(:,2)],'LineWidth',1);
xlabel('Number of trials');ylabel('Cronbachs alpha'); title('Cronbachs alpha on aggregate scores (converted to percentiles)')
set(gca,'XTick',1:2:length(minimum_trial_num:30),'XTickLabel',[minimum_trial_num:2:30])
xtickangle(0)
legend({'Percentiles','Raw aggregate scores'},'Location','best')
