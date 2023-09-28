data = readtable('../data/output/binary_responses/fear_binary_responses_only.csv');
trial_info = readtable('../data/output/binary_responses/fear_trial_info.csv');
trialN = size(data,2);
N = size(data,1);
emoNames = {'Angry','Fearful','Happy','Sad','Tender'};

difficulty_item = mean(data{:,1:trialN},'omitnan');
difficulty_p = mean(data{:,1:trialN}','omitnan');
rasch_mirt = readtable('../data/output/binary_responses/irt_models/rasch_mirt.csv');
pl2_mirt = readtable('../data/output/binary_responses/irt_models/2pl_mirt.csv');
pl3_mirt = readtable('../data/output/binary_responses/irt_models/3pl_mirt.csv');
participant_scores = readtable('../data/output/binary_responses/irt_models/participantScores.csv');
%% Evaluate Rasch model R Square-Accuracy
for k = 1:N
    for j = 1:trialN
        P_hat(k,j) = 1/(1+exp(-(participant_scores{k,1}-rasch_mirt{j,2})));
    end
end
y = table2array(data);
y_non_missing_idx = find(~isnan(y(:)));
y = y(y_non_missing_idx);
[B,~,~,~,STATS] = regress(y,[P_hat(y_non_missing_idx),ones(length(y),1)]);
P_binary = P_hat(y_non_missing_idx)>0.5;
accuracy_score = sum(P_binary==y)/length(y);
disp(['Accuracy Score: ' num2str(round(accuracy_score,2))]);
%% Accuracy plots
%find number of misclassifications per participant
for k=1:N
    y_p = data{k,:};
    y_idx = find(~isnan(y_p));
    binary_p = P_hat(k,y_idx)>0.5;
    accuracy_p(k) = sum(binary_p==y_p(y_idx))/length(binary_p);
end

%find number of misclassifications per item
for i=1:trialN
    y_p = data{:,i};
    y_idx = find(~isnan(y_p));
    binary_p = P_hat(y_idx,i)>0.5;
    accuracy_i(i) = sum(binary_p==y_p(y_idx))/length(binary_p);
end

%scatterplot classification accuracy-response accuracy
figure
subplot(1,2,1)
scatter(accuracy_i,difficulty_item,100,'filled');
xlabel('Classification Accuracy','FontSize',24);
ylabel('Response accuracy','FontSize',24);
title({['Item Classification Accuracy'], ['r=' num2str(round(corr(accuracy_i',...
    difficulty_item'),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on

subplot(1,2,2)
scatter(accuracy_p,difficulty_p,100,'filled');
xlabel('Classification Accuracy','FontSize',24);
ylabel('Response accuracy','FontSize',24);
title({['Participant Classification Accuracy'],...
    ['r=' num2str(round(corr(accuracy_p',difficulty_p'),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on

%scatterplot classification accuract-rasch model estimates
figure
subplot(1,2,1)
scatter(accuracy_i,rasch_mirt{:,2},100,'filled');
xlabel('Classification Accuracy','FontSize',24);
ylabel('Model difficulty estimates','FontSize',24);
title({['Item Classification Accuracy'], ['r=' num2str(round(corr(accuracy_i',...
    rasch_mirt{:,2}),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on

subplot(1,2,2)
scatter(accuracy_p,participant_scores{:,1},100,'filled');
xlabel('Classification Accuracy','FontSize',24);
ylabel('Model participant ability','FontSize',24);
title({['Participant Classification Accuracy'],...
    ['r=' num2str(round(corr(accuracy_p',participant_scores{:,1}),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on
%% Plots Rasch model
figure
subplot(1,2,1)
histogram(rasch_mirt{:,2});
ylabel('Count','FontSize',32);
xlabel('Θ (Item Difficulty)','FontSize',24);
title('Rasch Coefficients','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
xline(mean(rasch_mirt{:,2}),'-',{'Mean θ'},'LineWidth',3);

subplot(1,2,2)
histogram(participant_scores{:,1});
ylabel('Count','FontSize',32);
xlabel('Participant ability (ML)','FontSize',24);
title('Rasch ML estimates','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
xline(mean(participant_scores{:,1}),'-',{'Mean θ'},'LineWidth',3);

figure
scatter(rasch_mirt{:,2},difficulty_item,100,'filled')
rho = corr(rasch_mirt{:,2},difficulty_item');
xlabel('Θ',...
    'FontSize',24)
ylabel('Item Response Accuracy','FontSize',24)
title(['Item Response Accuracy and Rasch Coeffs: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on
%% Calculate True Test Scores
theta_step = 0.05;
theta_low = -5;
theta_high = 5;
theta_range = theta_low:theta_step:theta_high;
%rasch
for j = 1:trialN
    k=1;
    for th = theta_range
    p_theta(j,k) = 1/(1+exp(-rasch_mirt{j,1}*(th-rasch_mirt{j,2})));
    k = k+1;
    end
end

figure
plot(sum(p_theta)/trialN,'LineWidth',5)
ylabel('Test True Score','FontSize',32);
xlabel('Θ','FontSize',24);
set(gca,'XTick',1:20:length(theta_range),'XTickLabel',theta_low:theta_high)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on
%% Calculate Information and Standard error
%rasch
k=1;
for th = theta_range
    for j = 1:trialN
        p_correct = 1/(1+exp(-(th-rasch_mirt{j,2})));
        p_incorrect = 1-p_correct;
        information_test(j,k) = p_correct*p_incorrect;
        se_test(k) = 1./sqrt(sum(information_test(:,k)));
    end
    k = k+1;
end

figure
plot([sum(information_test)/trialN;se_test]','LineWidth',5)
ylabel('Estimate','FontSize',32);
xlabel('Θ','FontSize',24);
set(gca,'XTick',1:20:length(theta_range),'XTickLabel',theta_low:theta_high)
set(gca,'FontSize',32,'LineWidth',2)
title('Test Information & Standard Error','FontSize',36)
legend('Information','Standard Error','Location','best')
box on
grid on

%% Correlations of participants ability
difficulty_p = mean(data{:,1:trialN}','omitnan'); 
rho = corr([table2array(participant_scores),difficulty_p']);
figure
heatmap(rho,"FontSize",22)
ax = gca;
ax.XData = [participant_scores.Properties.VariableNames,{'Responce Accuracy'}]
ax.YData = [participant_scores.Properties.VariableNames,{'Responce Accuracy'}]
%% Plots 2PL model
figure
subplot(1,2,1)
histogram(pl2_mirt{:,1});
ylabel('Count','FontSize',32);
xlabel('Slope (Item Discrimination)','FontSize',24);
title('Alpha','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
xline(mean(pl2_mirt{:,1}),'-',{'Mean Slope'},'LineWidth',3);

subplot(1,2,2)
histogram(pl2_mirt{:,2});
ylabel('Count','FontSize',32);
xlabel('Intercept (Item Difficulty)','FontSize',24);
title('Beta','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
xline(mean(pl2_mirt{:,2}),'-',{'Mean Difficulty'},'LineWidth',3);


pl2_trials = find(pl2_mirt{:,1}>0.15)
figure
scatter(pl2_mirt{pl2_trials,2},difficulty_item(pl2_trials),100,'filled')
rho = corr(pl2_mirt{pl2_trials,2},difficulty_item(pl2_trials)');
xlabel('Θ',...
    'FontSize',24)
ylabel('Item Response Accuracy','FontSize',24)
title(['Item Response Accuracy and 2PL Coeffs: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on

%% Plots 3PL model
figure
subplot(1,2,1)
histogram(pl3_mirt{:,1});
ylabel('Count','FontSize',32);
xlabel('Slope (Item Discrimination)','FontSize',24);
title('3PL Model','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
xline(mean(pl3_mirt{:,1}),'-',{'Mean Slope'},'LineWidth',3);

subplot(1,2,2)
histogram(pl3_mirt{:,2});
ylabel('Count','FontSize',32);
xlabel('Intercept (Item Difficulty)','FontSize',24);
title('3PL Model','FontSize',28);
set(gca,'FontSize',32,'LineWidth',2)
xline(mean(pl3_mirt{:,2}),'-',{'Mean Diff.'},'LineWidth',3);

pl3_trials = find(pl2_mirt{:,1}>0.55)
figure
scatter(pl3_mirt{pl3_trials,2},difficulty_item(pl3_trials),100,'filled')
rho = corr(pl3_mirt{pl3_trials,2},difficulty_item(pl3_trials)');
xlabel('Θ',...
    'FontSize',24)
ylabel('Item Response Accuracy','FontSize',24)
title(['Item Response Accuracy and 3PL Coeffs: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on
