set(groot,'defaultAxesFontName','Arial')
data = readtable('data/output/binary_responses/binary_responses.csv');
trial_info = readtable('data/output/binary_responses/trial_info.csv');
trialN = size(data,2)-3;
emoNames = {'Anger','Fear','Happiness','Sadness','Tenderness'};
colors.blue = [0.1035, 0.4146, 0.6928];
colors.orange = [255, 127, 14] / 255;
difficulty_item = mean(data{:,1:trialN},'omitnan');
%% Boxplot of emotion item agreement
colorsfig1 = [
    0.1035, 0.4146, 0.6928;   % Blue
    0.9294, 0.6941, 0.1255;   % Dark yellow
    0.2431, 0.6824, 0.2431;   % Lighter green
    0.8392, 0.3294, 0.2196;   % Red-orange 
    0.5804, 0.4039, 0.7412;   % Purple
    0.25, 0.25, 0.25];
meanVals = grpstats(difficulty_item, trial_info{:,2}, 'median');
[~, sortIdx] = sort(meanVals, 'ascend');
all_agreements = [];
all_group_ids = [];
for i = sortIdx'
    idx = trial_info{:,2} == i;
    emo_vals = difficulty_item(idx);
    all_agreements = [all_agreements; emo_vals(:)];
    all_group_ids = [all_group_ids; repmat(i, length(emo_vals), 1)];
end
data_with_all = [all_agreements; all_agreements];
group_with_all = [all_group_ids; repmat(6, size(trial_info{:,2}))];
%figure('Position', [100, 100, 1000, 600])
figure
%t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
t = tiledlayout(1, 2);
nexttile(1)

hold on
% Plot scatter points
for i = sortIdx'
    this_vals = all_agreements(all_group_ids == sortIdx(i));
    jitter = (rand(size(this_vals)) - 0.5) * 0.3;
    scatter(this_vals, i + jitter, ...
        90, 'MarkerFaceColor', colorsfig1(i,:), ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
end
jitter = (rand(size(difficulty_item)) - 0.5) * 0.3; %scatter plot for all items
scatter(difficulty_item, 6 + jitter, ...
        90, 'MarkerFaceColor', 'k', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
% Plot horizontal boxplot
boxplot(data_with_all, categorical(group_with_all), ...
    'Orientation', 'horizontal', ...
    'Widths', 0.4, ...
    'Symbol', '','GroupOrder',cellstr([string(sortIdx);"6"]));
set(gca, 'FontSize', 22, 'LineWidth', 3)  
yticklabels([emoNames(sortIdx),{'All Items'}]);
xlabel('Item Agreement', 'FontSize', 24)
ylim([0.5, 6.5])
%title('Item Agreement by Emotion (Horizontal Boxplots)', 'FontSize', 24)
grid on
box on
text(-0.09,1.03,'a','Units','normalized','FontSize',30,'FontWeight','bold')
boxHandles     = flipud(findobj(gca, 'Tag', 'Box'));
whiskerHandles = flipud(findobj(gca, 'Tag', 'Whisker'));
medianHandles  = flipud(findobj(gca, 'Tag', 'Median'));
capHandles     = flipud(findobj(gca, 'Tag', 'Cap'));
allLines = findall(gca, 'Type', 'Line');

yVals = 1:6; 
yTol = 0.1;  

for k = 1:length(allLines)
    yData = get(allLines(k), 'YData');
    allLines(k).LineWidth = 2.5;
    if all(abs(yData - round(yData)) < yTol)
        allLines(k).LineStyle = '--';
    end
end
for i = 1:6
    ci = colorsfig1(i,:);
    if i <= length(boxHandles)
        set(boxHandles(i), 'Color', ci, 'LineWidth', 4)
    end
    if i <= length(medianHandles)
        set(medianHandles(i), 'Color', ci, 'LineWidth', 4)
    end
end
hold off
%% Create 1 emo pairs
X = table(categorical(trial_info{:,1}),categorical(trial_info{:,2}),difficulty_item');
emopair_cats = X.Var1;
g_mean_emopairs = groupsummary(X(:,[1,3]),'Var1','mean');
g_std_emopairs = groupsummary(X(:,[1,3]),'Var1','std');
g_mean_emopairs.std = g_std_emopairs{:,end};
g_mean_emopairs = sortrows(g_mean_emopairs,3,'descend');
%% Boxplot emo pairs
g_median_emopairs = groupsummary(X(:,[1,3]),'Var1','median');
g_median_emopairs = sortrows(g_median_emopairs,3,'descend');
% Generate distinct colors for each pair 
colors_pairs = turbo(20); 

% Verify colors are valid RGB triplets
assert(all(colors_pairs(:) >= 0 & colors_pairs(:) <= 1), 'Colors must be in [0,1]');

%figure('Position', [100, 100, 1400, 800])
nexttile(2)
hold on

for k = 1:20
    vals = difficulty_item(strcmpi(trial_info{:,1},string(g_median_emopairs{k,1})));
    jitter = (rand(size(vals)) - 0.5) * 0.3;
    scatter(vals, k + jitter, ...
        80, 'MarkerFaceColor', colors_pairs(k,:), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
end
boxplot(difficulty_item, trial_info{:,1}, ...
    'Orientation', 'horizontal', ...
    'Widths', 0.5, ...
    'Symbol', '','GroupOrder',cellstr(g_median_emopairs{:,1}));

boxHandles     = flipud(findobj(gca, 'Tag', 'Box'));
medianHandles  = flipud(findobj(gca, 'Tag', 'Median'));

allLines = findall(gca, 'Type', 'Line');

for k = 1:length(allLines)
    yData = get(allLines(k), 'YData');
    allLines(k).LineWidth = 2;
    if all(abs(yData - round(yData)) < yTol)
        allLines(k).LineStyle = '--';
    end
end

for k = 1:20
    ci = colors_pairs(k,:);
    if k <= length(boxHandles)
        set(boxHandles(k), 'Color', ci, 'LineWidth', 3)
    end
    if k <= length(medianHandles)
        set(medianHandles(k), 'Color', ci, 'LineWidth', 3)
    end
end

set(gca, 'FontSize', 22, 'LineWidth', 2,'YTick',1:20,'YTickLabel',g_median_emopairs{:,1})
ylim([0.5, 20 + 0.5])
text(-0.09,1.03,'b','Units','normalized','FontSize',30,'FontWeight','bold')
xlabel('Item Agreement', 'FontSize', 24)
ax = gca;              
ax.YAxis.FontSize = 18;
grid on
box on
hold off
%% Regressions
X_regress = X;
g_mean_emopairs = groupsummary(X(:,[1,3]),'Var1','mean');
%find reference category with minimum distance from mean
[m,i] = min(abs(g_mean_emopairs{:,3}-mean(difficulty_item)));
reference = g_mean_emopairs(i,1);
X_cats = categories(X_regress.Var1);
X_regress.Var1 = reordercats(X_regress.Var1,[X_cats(i);X_cats(1:i-1);X_cats(i+1:end)]);
model_emopairs = fitlm(X_regress,'Var3~Var1');
tbl = round(table2array(model_emopairs.Coefficients),2);
tbl = array2table(tbl,'VariableNames',model_emopairs.Coefficients.Properties.VariableNames,...
    'RowNames',[{'Intercept'};X_cats(1:i-1);X_cats(i+1:end)]);
writetable(tbl,'Table1,emopairs.csv','WriteRowNames',1)
model_targetemo = fitlm(X,'Var3~Var2');
%groupstats
g_mean_target = groupsummary(X(:,[2,3]),'Var2','mean');
g_std_target = groupsummary(X(:,[2,3]),'Var2','std');
g_mean_target.std = g_std_target{:,end};
g_mean_target = sortrows(g_mean_target,3,'descend');

%regression for demographic features
participant_scores = readtable('../data/output/binary_responses/irt_models/final_participantScores.csv');
X = table(data{:,trialN+1},data{:,trialN+2},...
    data{:,trialN+3},participant_scores{:,1});
%regression
model = fitlm(X,'Var4~Var1+Var2+Var3');
%% Rasch model plots
data = readtable('data/output/binary_responses/final_binary_responses.csv');
trial_info = readtable('data/output/binary_responses/final_trial_info.csv');
rasch_mirt = readtable('data/output/binary_responses/irt_models/final_rasch_mirt.csv');
participant_scores = readtable('data/output/binary_responses/irt_models/final_participantScores.csv');

N = size(data,1);
guessing = 0.5;
trialN = size(data,2);
difficulty_item = mean(data{:,1:trialN},'omitnan');
difficulty_p = mean(data{:,1:trialN}','omitnan');
%% Calculate TSS and Information
%Calculate True Test Scores
theta_step = 0.05;
theta_low = -5;
theta_high = 5;
theta_range = round(theta_low:theta_step:theta_high,2);
%rasch
for j = 1:trialN
    k=1;
    for th = theta_range
        p_theta(j,k) = guessing+(1-guessing)*(1/(1+exp(-rasch_mirt{j,1}*(th-rasch_mirt{j,2}))));
        k = k+1;
    end
end

%TSS and Information
k=1;
for th = theta_range
    for j = 1:trialN
        p_correct(j,k) = guessing+(1-guessing)*(1/(1+exp(-(th-rasch_mirt{j,2}))));
        p_incorrect(j,k) = 1-p_correct(j,k);
        product = p_correct(j,k)*p_incorrect(j,k);
        if k>1
            j_der1 = [p_correct(j,k)-p_correct(j,k-1)]./theta_step;
            information_test(j,k) = j_der1^2/product;
        else
            information_test(:,1) =  NaN;
        end
    end
    se_test(k) = 1./sqrt(sum(information_test(:,k)));
    k = k+1;
end
information_test(:,1) =  information_test(:,2);
se_test(1) = 1./sqrt(sum(information_test(:,1)));
tss = sum(p_theta)/trialN;
tss_sc = tss(theta_range<2.4 & theta_range>-2.4);%scatter plot axis for figure 2
x_sc = theta_range(theta_range<2.4 & theta_range>-2.4);
%% Plots IRT model
%colors.blue = [44,101,220] / 255;
colors.blue = brewermap(10,'Blues');
colors.blue = colors.blue(8,:);
colors.orange = [255, 127, 14] / 255;
figure
ax1 = subplot(2,2,1)
hold on
histogram(rasch_mirt{:,2},'Normalization','pdf','FaceColor',colors.blue,'LineWidth',1.5,'BinEdges',-5.5:5.5);
[f, x] = ksdensity(rasch_mirt{:,2});
xlim([-5,5])
set(gca,'XTick',-5:5)
set(gca,'FontSize',28,'LineWidth',2)
plot(x, f, 'LineWidth', 4, 'Color', colors.orange)
ylabel('PDE');
xlabel('Model difficulty estimates (b)');
title('Item Difficulty');
text(-0.2,1,'a','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
hold off

ax2 = subplot(2,2,2)
hold on
histogram(participant_scores{:,1},'Normalization','pdf','FaceColor',colors.blue,'LineWidth',1.5,'BinEdges',-5.5:5.5);
[f, x] = ksdensity(participant_scores{:,1});
xlim([-5,5])
set(gca,'XTick',-5:5)
set(gca,'FontSize',28,'LineWidth',2)
plot(x, f, 'LineWidth', 4, 'Color', colors.orange )
ylabel('PDE');
xlabel('Model ability estimates (θ)');
title('Participant Ability');
text(-0.2,1,'b','Units','normalized','FontSize',30,'FontWeight','bold')
%xline(mean(participant_scores{:,1}),'-',{'Mean θ'},'LineWidth',3);
box on
grid on
hold off
%linkaxes([ax1, ax2], 'y');

ax3 = subplot(2,2,3)
plot(round(sum(p_theta)/trialN,4),'LineWidth',4,'Color',colors.blue)
set(gca,'FontSize',28,'LineWidth',2)
title('Test True Score');
ylabel('Correct Responses')
xlabel('Ability Level (θ)');
set(gca,'XTick',1:20:length(theta_range),'XTickLabel',theta_low:theta_high)
xlim([1,length(theta_range)])
text(-0.18,1,'c','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on

ax4 = subplot(2,2,4)
hold on
plot(round(rescale(sum(information_test)/trialN),4)','LineWidth',4,'Color',colors.blue)
plot(round(rescale(se_test),4)','LineWidth',4,'Color',colors.orange)
set(gca,'FontSize',28,'LineWidth',2)
ylabel('Estimate');
xlabel('Ability Level (θ)');
set(gca,'XTick',1:20:length(theta_range),'XTickLabel',theta_low:theta_high)
title('Test Information & Standard Error')
legend('Information','Standard Error','Location','best')
xlim([1,length(theta_range)])
text(-0.18,1,'d','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
hold off
%linkaxes([ax3, ax4], 'y');
set(ax1,'Position',[0.13,0.6,0.3347,0.3412])
set(ax2,'Position',[0.5703,0.6,0.3347,0.3412])
set(ax3,'Position',[0.13,0.0950,0.3347,0.3412])
set(ax4,'Position',[0.5703,0.0950,0.3347,0.3412])
%% Rasch model Figure 3
% Evaluate Rasch model R Square-Accuracy
for k = 1:N
    for j = 1:trialN
        P_hat(k,j) = guessing+(1-guessing)*(1/(1+exp(-(participant_scores{k,1}-rasch_mirt{j,2}))));
    end
end

y = table2array(data);
y_non_missing_idx = find(~isnan(y(:)));
y = y(y_non_missing_idx);
bs = mean((P_hat(y_non_missing_idx)-y).^2); %brier score formula
P_binary = P_hat(y_non_missing_idx)>0.75;
accuracy_score = sum(P_binary==y)/length(y);
%% Figure 3 with subplots
scatter_color = [0.2228, 0.5277, 0.7542];
orange_trans = [1.0000 0.4980 0.0549 0.9];
orange_ci = [1.0000 0.4980 0.0549 0.4];
%find number of misclassifications per participant
for k=1:N
    y_p = data{k,:};
    y_idx = find(~isnan(y_p));
    binary_p = P_hat(k,y_idx)>0.75;
    accuracy_p(k) = sum(binary_p==y_p(y_idx))/length(binary_p);
end

%find number of misclassifications per item
for i=1:trialN
    y_p = data{:,i};
    y_idx = find(~isnan(y_p));
    binary_p = P_hat(y_idx,i)>0.75;
    accuracy_i(i) = sum(binary_p==y_p(y_idx))/length(binary_p);
end

figure
subplot(2,2,[1 2])
hold on
scatter(rasch_mirt{:,2},difficulty_item,80,scatter_color,'filled')
[x_unique, ~, idx] = unique(rasch_mirt{:,2});
y_mean = accumarray(idx, difficulty_item, [], @mean);
y_smooth = smooth(x_unique, y_mean, 0.3, 'loess');
x_fit = linspace(min(x_unique), max(x_unique), 200)';
y_fit = interp1(x_unique, y_smooth, x_fit, 'pchip');
% Adaptive confidence intervals
resid = y_mean - y_smooth;  % residuals
bandwidth = 1;              % controls smoothness of CI (bigger = smoother)
local_sigma = zeros(size(x_fit));

for i = 1:length(x_fit)
    % Gaussian weights for local variance estimation
    w = exp(-0.5 * ((x_unique - x_fit(i)) / bandwidth).^2);
    w = w / sum(w);
    local_sigma(i) = sqrt(sum(w .* resid.^2));
end

% 95% CI
ci = 1.96 .* local_sigma;
y_upper = y_fit + ci;
y_lower = y_fit - ci;
fill([x_fit; flipud(x_fit)], [y_upper; flipud(y_lower)], ...
    [0.2 0.6 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_fit,y_fit, 'Color',orange_trans, 'LineWidth', 4);
%set(gca,'XTick',-5:5)
set(gca,'FontSize',26,'LineWidth',2)
%rho = corr(rasch_mirt{:,2},difficulty_item');
xlabel('Item Difficulty (Model)')
ylabel('Item Agreement')
title('Item Agreement and Item Difficulty')
text(-0.1,1,'a','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
hold off

ax3 = subplot(2,2,3)
hold on
s=scatter(rasch_mirt{:,2},accuracy_i,80,scatter_color,'filled');
[x_unique, ~, idx] = unique(rasch_mirt{:,2});
y_mean = accumarray(idx, accuracy_i, [], @mean);
y_smooth = smooth(x_unique, y_mean, 0.3, 'loess');
x_fit = linspace(min(x_unique), max(x_unique), 200)';
y_fit = interp1(x_unique, y_smooth, x_fit, 'pchip');
%Adaptive confidence intervals
resid = y_mean - y_smooth;  % residuals
bandwidth = 1;              % controls smoothness of CI (bigger = smoother)
local_sigma = zeros(size(x_fit));

for i = 1:length(x_fit)
    % Gaussian weights for local variance estimation
    w = exp(-0.5 * ((x_unique - x_fit(i)) / bandwidth).^2);
    w = w / sum(w);
    local_sigma(i) = sqrt(sum(w .* resid.^2));
end

% 95% CI
ci = 1.96 .* local_sigma;
y_upper = y_fit + ci;
y_lower = y_fit - ci;
fill([x_fit; flipud(x_fit)], [y_upper; flipud(y_lower)], ...
    [0.2 0.6 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_fit,y_fit, 'Color',orange_trans, 'LineWidth', 4);
yticks([0:0.2:1]);
xticks([-5:5]);
set(gca,'FontSize',26,'LineWidth',2)
ylabel('Accuracy');
xlabel('Item difficulty (Model)');
title('Accuracy and Item Difficulty')
%title({['Item Classification Accuracy'], ['r=' num2str(round(corr(accuracy_i',...
%    rasch_mirt{:,2}),2))]},'FontSize',28)
text(-0.25,1,'b','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
hold off

ax4 = subplot(2,2,4)
hold on
scatter(participant_scores{:,1},accuracy_p,80,scatter_color,'filled');
y_smooth = smooth(participant_scores{:,1}, accuracy_p, 0.5, 'loess');
x_fit = linspace(min(participant_scores{:,1}), max(participant_scores{:,1}), 200)';
y_fit = interp1(participant_scores{:,1}, y_smooth, x_fit, 'pchip');
%Adaptive confidence intervals
bandwidth = 1;  % controls smoothness of CI (larger = smoother, wider bands)
resid = accuracy_p' - y_smooth;
local_sigma = zeros(size(x_fit));

%% HERE
for i = 1:length(x_fit)
    % Gaussian weights based on distance
    w = exp(-0.5 * ((participant_scores{:,1} - x_fit(i)) / bandwidth).^2);
    w = w / sum(w);  % normalize
    % Weighted std of residuals
    local_sigma(i) = sqrt(sum(w .* resid.^2));
end

% 95% pointwise CI
ci = 1.96 .* local_sigma;
y_upper = y_fit + ci;
y_lower = y_fit - ci;
fill([x_fit; flipud(x_fit)], [y_upper; flipud(y_lower)], ...
    [0.2 0.6 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(x_fit,y_fit, 'Color',orange_trans, 'LineWidth', 4);
yticks([0:0.2:1]);
xticks([-5:5]);
xlim([-2.7,2.4])
set(gca,'FontSize',26,'LineWidth',2)
%plot(x_sc,tss_sc)
ylabel('Accuracy');
xlabel('Participant ability (Model)');
title('Accuracy and Participant Ability')
%title({['Participant Classification Accuracy'],...
%    ['r=' num2str(round(corr(accuracy_p',participant_scores{:,1}),2))]},'FontSize',28)
text(-0.15,1,'c','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
hold off
set(ax3,'Position',[0.13,0.0950,0.3347,0.3412])
set(ax4,'Position',[0.5703,0.0950,0.3347,0.3412])