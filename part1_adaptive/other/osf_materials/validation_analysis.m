%% Validation analysis 
data = readtable('data/validation/study3Data.csv');
theta = table2array(readtable('data/validation/study3Theta.csv'));
difficulty = table2array(readtable('data/validation/difficulty.csv'));
%% Calculate Information and SE
guessing = 0.5;
theta_step = 0.02;
for k = 1:size(theta,1)
    for j = 1:20
        count = 1;
        temp_information_test = 0;
        while count<=j
        p_correct = guessing+(1-guessing)*(1/(1+exp(-(theta(k,j)-difficulty(k,count)))));
        p_incorrect = 1-p_correct;
        product = p_correct*p_incorrect;
        temp_information_test = temp_information_test+(product/0.5^2);
        count = count+1;
        end
        information_test(k,j) = temp_information_test;
        se(k,j) = 1/(sqrt(information_test(k,j)));
        rel(k,j) = 1 - se(k,j)^2;
    end
end
colors.blue = [0.1035, 0.4146, 0.6928];
colors.orange = [1 0.4980 0.0549,0.85];
n = size(se,1);
alpha = 0.05;
z = norminv(1 - alpha/2);
for i = 1:20
    mean_se(i) = nanmean(se(:,i));
    mean_rel(i) = nanmean(rel(:,i));
    std_se(i) = nanstd(se(:,i));
    %calculate assymetric CI for se
    se_sorted_data = sort(se(:,i));
    alpha = 0.05;
    lower_percentile = alpha / 2;
    upper_percentile = 1 - (alpha / 2);
    se_lower_bound(i) = prctile(se_sorted_data, lower_percentile * 100);
    se_upper_bound(i) = prctile(se_sorted_data, upper_percentile * 100);
    %calculate assymetric CI for reliability
    rel_sorted_data = sort(rel(:,i));
    alpha = 0.05;
    rel_lower_bound(i) = prctile(rel_sorted_data, lower_percentile * 100);
    rel_upper_bound(i) = prctile(rel_sorted_data, upper_percentile * 100);
end
%% Figure 7 ability histogram and SE
for i = 1:20
    %se
    mean_se(i) = nanmean(se(:,i));
    se_sd = nanstd(se(:,i));
    se_lower_bound(i) = mean_se(i) - z * (se_sd / sqrt(n));
    se_upper_bound(i) = mean_se(i) + z * (se_sd / sqrt(n));

    %reliability
    mean_rel(i) = nanmean(rel(:,i));
    rel_sd = nanstd(rel(:,i));
    rel_lower_bound(i) = mean_rel(i) - z * (rel_sd / sqrt(n));
    rel_upper_bound(i) = mean_rel(i) + z * (rel_sd / sqrt(n));

end

figure
subplot(1,2,1)
hold on
histogram(theta(:,end),[-5:5],'Normalization','pdf','FaceColor',colors.blue,'LineWidth',1.5);
[f, x] = ksdensity(theta(:,end));
plot(x, f, 'LineWidth', 4, 'Color', colors.orange)
set(gca,'FontSize',28,'LineWidth',2)
[f, x] = ksdensity(theta(:,end));
%set(gca,'XTick',-5:5)
ylabel('PDE');
xlabel('SMART Score (θ)');
set(gca,'XTick',-5:5)
xlim([-5,5])
mean_t = mean(theta(:,end));
title('SMART Score');
text(-0.2,1,'a','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
hold off

subplot(1,2,2)
hold on;
plot(1:length(mean_se),mean_se,'LineWidth', 5, 'Color', colors.blue)
plot(1:length(mean_rel),mean_rel,'LineWidth', 5, 'Color', colors.orange)
fill([1:length(mean_se), fliplr(1:length(mean_se))], [se_upper_bound, fliplr(se_lower_bound)], [0.2 0.6 0.8], ...
    'FaceAlpha', 0.5,'EdgeColor','none'); 
set(gca,'FontSize',28,'LineWidth',2)
fill([1:length(mean_rel), fliplr(1:length(mean_rel))], [rel_upper_bound, fliplr(rel_lower_bound)], [1 0.65 0.3], ...
    'FaceAlpha', 0.5,'EdgeColor','none'); 
xlim([2,20])
set(gca,'XTick',0:2:20)
ylabel('Estimate')
%ylim([-1.5,1.5])
%set(gca,'YTick',-1.5:0.25:1.5)
xlabel('Test Length')
title('SEM & Conditional Reliability')
text(-0.15,1,'b','Units','normalized','FontSize',30,'FontWeight','bold')
legend({'Standard Error', 'Reliability'},'FontSize',24)
box on 
grid on
hold off
%% Regression
predictors = [5,6,7,15,16,18,19];

predictors_names = {'Tas_Total','MSI_MTraining','mProms_total','score_facial',...
    'score_vocal','Age','Gender'};
formula = ['Theta ~ Tas_Total + MSI_MTraining + mProms_total + score_facial',...
    +' + score_vocal + Age + Sex_bin'];
%zscore
data_reg = data;
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
data_reg{:,[1 predictors]} = zscor_xnan(data{:,[1 predictors]});
mdl = fitlm(data_reg(:,[1 predictors]),formula);
disp('Linear regression formula')
disp(formula)
disp(['R Square: ',num2str(round(mdl.Rsquared.Ordinary,2))])
mdl.Coefficients(:,[1,4])
%% Correlation matrix
vars = [1,5,6,7,15,16,18];
rho_main = corr(table2array(data(:,vars)),'rows','pairwise');
figure
imagesc(rho_main);colorbar();
set(gca,'YTick',1:length(vars),'YTickLabel',data.Properties.VariableNames(vars))
set(gca,'XTick',1:length(vars),'XTickLabel',data.Properties.VariableNames(vars))

rho = corr(table2array(data(:,1:16)),'rows','pairwise');
figure
imagesc(rho);colorbar();
set(gca,'YTick',1:19,'YTickLabel',data.Properties.VariableNames)
set(gca,'XTick',1:19,'XTickLabel',data.Properties.VariableNames)
%% Plots
%Alexithymia plot
alex_ratio(1) = sum(data.Tas_Total<=51)/height(data);
alex_ratio(3) = sum(data.Tas_Total>=61)/height(data);
alex_ratio(2) = sum(data.Tas_Total<61 & data.Tas_Total>51)/height(data);
figure
bar(alex_ratio)
set(gca,'XTick',1:3,'XTickLabel',{'No Alexithymia','Possible Alexithymia',...
    'Alexithymia'},'FontSize',18)
%ylabel('Percentage','FontSize',20);
title('TAS results')

%microproms
figure
subplot(2,2,1)
histogram(data.mProms_total,8)
set(gca,'FontSize',18)
xlabel('MicroProms')
xline(11.96,'-',{'Original study mean'},'LineWidth',3)
title('MicroProms')
%msi
subplot(2,2,2)
histogram(data.MSI_MTraining,8)
set(gca,'FontSize',18)
xlabel('Musical Training')
xline(26.5,'-',{'Original study mean'},'LineWidth',3)
title('MSI Scores')

%ERI
subplot(2,2,3)
histogram(data.score_facial,8)
set(gca,'FontSize',18)
xlabel('ERI facial')
xline(71,'-',{'Original study mean'},'LineWidth',3)
title('ERI facial Scores')

subplot(2,2,4)
histogram(data.score_vocal,8)
set(gca,'FontSize',18)
xlabel('ERI vocal')
xline(71,'-',{'Original study mean'},'LineWidth',3)
title('ERI vocal Scores')

%Participants age
figure
histogram(data.Age,8);
set(gca,'FontSize',18)
xlabel('Age');
title('Participant Age');
%% Make correlation table
vars_to_compare = {[6,15,16,18],[7:14],[2:5]};
significance = {[0.21,0.26,0.31],[0.189,0.232,0.279]};
for k = 1:3
    n = 1;
    for v = vars_to_compare{k}
        N_tbl{k}(n) = sum(~isnan(data{:,v}));
        mean_tbl{k}(n) = nanmean(data{:,v});
        std_tbl{k}(n) = nanstd(data{:,v});
        if k ==3
            [rho_tbl{k}(n),pval{k}(n)] = corr(theta(:,end),data{:,v},'rows','pairwise','tail','left');
        else
            [rho_tbl{k}(n),pval{k}(n)] = corr(theta(:,end),data{:,v},'rows','pairwise','tail','right');
        end
        n = n+1;
    end
end
corr_tbl.names = [data.Properties.VariableNames(vars_to_compare{1}),data.Properties.VariableNames(vars_to_compare{2}),...
    data.Properties.VariableNames(vars_to_compare{3})]'
corr_tbl.N = [N_tbl{1},N_tbl{2},N_tbl{3}]';
corr_tbl.M = round([mean_tbl{1},mean_tbl{2},mean_tbl{3}]',2);
corr_tbl.SD = round([std_tbl{1},std_tbl{2},std_tbl{3}]',2);
corr_tbl.Rho = round([rho_tbl{1},rho_tbl{2},rho_tbl{3}]',2);
corr_tbl.p = round([pval{1},pval{2},pval{3}]',3);
corr_tbl=struct2table(corr_tbl)
%writetable(corr_tbl,'../paper/corr_table.csv')
%% Percentage of data within -+ 1SD
upper_thx = nanmean(theta(:,end)) + nanstd(theta(:,end));
low_thx = nanmean(theta(:,end)) - nanstd(theta(:,end));
percentage = sum(theta(:,end)<upper_thx & theta(:,end)>low_thx)/size(theta,1);