%publish('pilot_irt_conversion.m','format','pdf','showCode',false);
%cd Documents/projects/github/MMBB_emotion_demo/part1_adaptive/other/analysis/
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
data = readtable('pilot_data.xlsx');
lvl_achieved_col = 121;
difficulty_col = 61;
level_ended = cell(1,2); difficulty_sum = cell(1,2); lvl_achieved = cell(1,2);
minimum_trial_num = 12;
trials = readtable('../data/output/trials.csv');
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
xlabel('Trials');ylabel('TStatistic');title('Trial difficulty per emotion');
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

%%
figure
scatter(1-rescale(trials.Distance(4:end)),occurences(4:end),80,'filled')
xlabel('Item difficulty','FontSize',32); 
ylabel('Item instances','FontSize',32);
set(gca,'FontSize',32,'LineWidth',2)
grid on
box on
ylim([0 18])
%title('Number of instances per trial')
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
x = 1-rescale(trials.Distance(~isnan(trial_response_ratio)));
y = trial_response_ratio(~isnan(trial_response_ratio))';
%fit linear model
[b,~,~,~,stats] = regress(y,[x,ones(size(x,1),1)]);
p = polyfit(x,y,2);
p3 = polyfit(x,y,3);
y_polyval = polyval(p,x);
y3_polyval = polyval(p3,x);

USS = sum((y-y_polyval).^2);
TSS = sum((y-mean(y)).^2);
ESS =  TSS - USS;
r_square(1) = ESS/TSS;

USS = sum((y-y3_polyval).^2);
TSS = sum((y-mean(y)).^2);
ESS =  TSS - USS;
r_square(2) = ESS/TSS;
x_plot = [0:0.01:1];

figure
subplot(1,3,1)
hold on
scatter(1-rescale(trials.Distance),trial_response_ratio,80,'filled');
xlabel('Difficulty','FontSize',20);ylabel('Correct response ratio','FontSize',20);
title({'Linear model',['Model: y = ' , num2str(round(b(1),2)),...
    'x + ', num2str(round(b(2),2))],['R^2 = ' num2str(round(stats(1),2))]},'FontSize',20);
plot([0:1],b(1)*[0:1]+b(2),'LineWidth',3)
axis([0.2 1 0 1.2])
set(gca,'FontSize',20,'LineWidth',2)
box on
grid on
hold off

subplot(1,3,2)
hold on
scatter(1-rescale(trials.Distance),trial_response_ratio,80,'filled');
xlabel('Difficulty','FontSize',20);ylabel('Correct response ratio','FontSize',20);
title({'Quadratic',['Model: y = ' , num2str(round(p(1),2)),...
    'x^2 + ', num2str(round(p(2),2)), 'x + ', num2str(round(p(3),2))],...
    ['R^2 = ' num2str(round(r_square(1),2))]},'FontSize',20);
plot(x_plot,p(1)*(x_plot.^2)+p(2)*x_plot+p(3),'LineWidth',3)
axis([0.2 1 0 1.2])
set(gca,'FontSize',20,'LineWidth',2)
box on
grid on
hold off

subplot(1,3,3)
hold on
scatter(1-rescale(trials.Distance),trial_response_ratio,80,'filled');
xlabel('Difficulty','FontSize',20);ylabel('Correct response ratio','FontSize',20);
title({'Cubic',['Model: y = ' , num2str(round(p3(1),2)),'x^3 + ',num2str(round(p3(2),2)),...
    'x^2 + ', num2str(round(p3(3),2)), 'x + ', num2str(round(p3(4),2))],...
    ['R^2 = ' num2str(round(r_square(2),2))]},'FontSize',20);
plot(x_plot,p3(1)*(x_plot.^3)+p3(2)*(x_plot.^2)+p3(3)*x_plot+p3(4),'LineWidth',3)
axis([0.2 1 0 1.2])
set(gca,'FontSize',20,'LineWidth',2)
box on
grid on
hold off

y_cubic_pred = p3(1)*(x_plot.^3)+p3(2)*(x_plot.^2)+p3(3)*x_plot+p3(4);