data = readtable('../data/output/binary_responses/fear_binary_responses_only.csv');
trial_info = readtable('../data/output/binary_responses/fear_trial_info.csv');
trialN = size(data,2);
emoNames = {'Angry','Fearful','Happy','Sad','Tender'};
%% Analysis of emotion types
difficulty = mean(data{:,1:trialN},'omitnan');
X = table(trial_info{:,1},categorical(trial_info{:,2}),difficulty');
%regressions
model_emopairs = fitlm(X,'Var3~Var1');
model_targetemo = fitlm(X,'Var3~Var2');
%groupstats
g_mean_target = groupsummary(X(:,[2,3]),'Var2','mean');
g_std_target = groupsummary(X(:,[2,3]),'Var2','std');
g_mean_target.std = g_std_target{:,end};
g_mean_target = sortrows(g_mean_target,3,'descend');
g_mean_emopairs = groupsummary(X(:,[1,3]),'Var1','mean');
g_std_emopairs = groupsummary(X(:,[1,3]),'Var1','std');
g_mean_emopairs.std = g_std_emopairs{:,end};
g_mean_emopairs = sortrows(g_mean_emopairs,3,'descend');
%indices of emo pairs anger-fear, fear-anger
combos = [1,5;2,9;3,13;4,17;6,10;7,14;8,18;11,15;12,19;16,20];
% Plots
figure
%positions = [1:height(g_mean_emopairs)];
errorbar(g_mean_emopairs{:,3},2:height(g_mean_emopairs)+1,g_mean_emopairs{:,4}*0.5,'horizontal','*','LineWidth',5,...
    'MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gca,'linewidth',3) 
set(gca,'ytick',2:height(g_mean_emopairs)+1,'yticklabels',...
    g_mean_emopairs{:,1},'Fontsize',30)
ylim([1 height(g_mean_emopairs)+2])
xlabel('Response Accuracy','FontSize',36)
title('Response Accuracy across EmoPairs','FontSize',36)

figure
%positions = [1:height(g_mean_emopairs)];
errorbar(g_mean_target{:,3},2:height(g_mean_target)+1,g_mean_target{:,4}*0.5,'horizontal','*','LineWidth',5,...
    'MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gca,'linewidth',3) 
set(gca,'ytick',2:height(g_mean_target)+1,'yticklabels',...
    emoNames(g_mean_target{:,1}),'Fontsize',30)
ylim([1 height(g_mean_target)+2])
xlabel('Response Accuracy','FontSize',36)
title('Response Accuracy across Emotions','FontSize',36)

figure
subplot(2,1,1)
plot(sort(difficulty),'LineWidth',5)
ylabel('Response accuracy','FontSize',32);
xlabel('Items','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(difficulty)])
box on
grid on

subplot(2,1,2)
hold on
for i = 1:5
    plot(sort(difficulty(find(trial_info{:,2}==i))),'LineWidth',5)
end
ylabel('Response accuracy','FontSize',32);
xlabel('Items','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 sum(trial_info{:,2}==mode(trial_info{:,2}))])
box on
grid on
legend(emoNames,'Location','best')
hold off

%% Analysis of demographic features
difficulty = mean(data{:,1:trialN}','omitnan');
X = table(data{:,trialN+1},data{:,trialN+2},...
    data{:,trialN+3},difficulty');
%regression
model = fitlm(X,'Var4~Var1+Var2+Var3');
%plots participant response accuracy per country
figure
subplot(2,1,1)
plot(sort(difficulty),'LineWidth',5)
ylabel('Response accuracy','FontSize',32);
xlabel('Participants','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(difficulty)])
box on
grid on
title('Response Accuracy per Participant','FontSize',36)

subplot(2,1,2)
hold on
plot(sort(difficulty(find(strcmpi(data{:,trialN+1},'fi')))),'LineWidth',5)
plot(sort(difficulty(find(strcmpi(data{:,trialN+1},'spa')))),'LineWidth',5)
hold off
ylabel('Response accuracy','FontSize',32);
xlabel('Participants','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 130])
box on
grid on
title('Response Accuracy Finnish-Spanish','FontSize',36)
legend({'Finnish','Spanish'},'Location','best')

%age plot
figure
boxplot(sort(data.Age))
ylabel('Age','FontSize',32);
xlabel('Participants','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(difficulty)])
box on
grid on
title('Age (Finnish participants)','FontSize',36)

figure
[rho,pval] = corr(difficulty',data.Age,'rows','pairwise');
scatter(difficulty',data.Age,150,'filled')
xlabel('Response Accuracy',...
    'FontSize',24)
ylabel('Age','FontSize',24)
title(['Scatterplot of Age and Response Accuracy: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on

[h,p,~,stats] = ttest2(difficulty(find(strcmpi(data.Gender,'Male'))),...
    difficulty(find(strcmpi(data.Gender,'Female'))))
figure
h = boxplot(difficulty,data.Gender,'OutlierSize',16);
set(h,{'linew'},{3})
ylabel('Response Accuracy','FontSize',32);
xlabel('Gender','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on
title(['Gender and Response Accuracy: tstat=',num2str(round(stats.tstat,2)),...
    ' pval=', num2str(round(p,3))],'FontSize',36)
hold off

%%Compare theta and item values
figure
subplot(1,2,1)
histogram(mean(data{:,1:trialN},'omitnan'));
ylabel('Count','FontSize',32);
xlabel('Response Accuracy','FontSize',24);
title('Item Difficulty','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)

subplot(1,2,2)
histogram(mean(data{:,1:trialN}','omitnan'));
ylabel('Count','FontSize',32);
xlabel('Response Accuracy','FontSize',24);
title('Participant Ability','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)