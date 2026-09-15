% HARDCODED MICROPROMS-CHECK DATA FILE
%publish('validation_analysis.m','format','pdf','showCode',false);
warning('off')
jspsych_dir = '../validation_study/data/jatos_results/';
jatos_metadata = '../validation_study/data/jspsych_metadata.json';
pilot_participants = [23:31,33,34];
duplicate_participants = {'60fe9a0daa398bbc0d610dc0'};
pilot_ids = {'662fc8f6aa3d3b206672e28d','5c181b8d905b250001be42ca','66a27a561b8ed92d82f3b97a',...
    '5ea1a2e7939362055f0c325a','5dc9c2c0b9fad36ddca37632','663a3ee80ab62f61353a10f0','562a00dac8ffc20012513fbe',...
    '55d22025cc2b18000c0b9d9c','6715521230e88bbd98ef83ba','62d6d87d91cf24a3c1c0f10c','5d49d17b3dad1f0001e2aba1'};
%% Load jspsych data
%53 = MSI1,2
%54 = MSI3:7
%55 = TAS1:10
%56 = TAS11:20
%load jatos metadata
jatos_meta = fileread(jatos_metadata);
meta = jsondecode(jatos_meta);
for k = 1:length(meta.data.studyResults)
    participant_starttime(k,:) = [meta.data.studyResults{k}.id,...
        meta.data.studyResults{k}.startDate];
end

filenames = dir(jspsych_dir);
filenames = filenames(4:end); %HARDCODED FOR PILOT
for i = 1:length(filenames)
    jspsych_table{i} = table();
    filename_number = strsplit(filenames(i).name,'_');
    filename_number = filename_number{3};
    filename = [jspsych_dir, filenames(i).name, '/comp-result_',filename_number,'/data.txt'];
    r = fileread(filename);
    if any(pilot_participants==str2num(filename_number))
        r = r(1:length(r)/2); %HARDCODED FOR PILOT
    end
    jspsych_data = jsondecode(r);
    counter = 0;
    if strcmpi(jspsych_data.trials{4}.response.audioTest,'orange')
       trial_idx = [13,53,54,55,56];
    else
       trial_idx = [13,53,54,55,56]+1; 
    end
    for t=1:length(jspsych_data.trials)
        if strcmpi(jspsych_data.trials{t}.trial_type,'audio-button-response') & t>=trial_idx(1)
            counter = counter+1;
            jspsych_table{i}{counter,'rt'} = jspsych_data.trials{t}.rt;
            jspsych_table{i}{counter,'trial_index'} = jspsych_data.trials{t}.trial_index;
            jspsych_table{i}{counter,'response'} = jspsych_data.trials{t}.response;
            jspsych_table{i}{counter,'Correct_emotion'} = string(jspsych_data.trials{t}.Correct_emotion);
            jspsych_table{i}{counter,'Wrong_emotion'} = string(jspsych_data.trials{t}.Wrong_emotion);
            jspsych_table{i}{counter,'Trial'} = jspsych_data.trials{t}.Trial;
            jspsych_table{i}{counter,'Track'} = jspsych_data.trials{t}.Track;
            jspsych_table{i}{counter,'Theta'} = jspsych_data.trials{t}.Theta;
            jspsych_table{i}{counter,'Difficulty'} = jspsych_data.trials{t}.Difficulty;
            jspsych_table{i}{counter,'Correct_response'} = jspsych_data.trials{t}.Correct_response;
            jspsych_table{i}{counter,'Optimizer'} = jspsych_data.trials{t}.Optimizer;
            if jspsych_data.trials{t}.Optimizer==2
                jspsych_table{i}{counter,'Iterations'} = jspsych_data.trials{t}.Iterations;
            end
            jspsych_table{i}{counter,'Epoch'} = jspsych_data.trials{t}.Epoch;
        end
    end
    if any(pilot_participants==str2num(filename_number))
        userID{i} = pilot_ids{find(pilot_participants==str2num(filename_number))};
    else
        try
            userID{i} = jspsych_data.trials{1}.userID;
        catch
            userID{i} = '00';
        end
    end
    % add study start time
    starttime(i) = participant_starttime(find(participant_starttime(:,1)==str2num(filename_number)),2);
end
for i = 1:length(userID) %drop invalid ids
    valid_id(i)=length(userID{i})>10;
end
for i = 1:length(duplicate_participants)
if any(strcmpi(userID,duplicate_participants{i}))
   userID{find(strcmpi(userID,duplicate_participants{i}),1)}=[duplicate_participants{i},'_dble']; 
end
end

%create matrix of thetas
for k = 1:(length(jspsych_table))
    theta(k,:) = jspsych_table{k}.Theta;
    difficulty(k,:) = jspsych_table{k}.Difficulty;
    response(k,:) = jspsych_table{k}.response;
    trial(k,:) = jspsych_table{k}.Trial;
end

%Drop non prolific cases
theta = theta(valid_id,:);
userID = userID(valid_id);
trial_difficulty = [];
mean_correct = [];
N = [];
counter = 0;
for i = unique(trial(:))'
    counter = counter+1;
    idx = find(i==trial);
    trial_difficulty(counter) = difficulty(idx(1));
    mean_correct(counter) = nanmean(response(idx));
    N(counter) = length(idx);
end
trial_difficulty = trial_difficulty(find(N>20));
mean_correct = mean_correct(find(N>20));
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
tiledlayout(1,2,'TileSpacing','loose','Padding','loose')
nexttile;
hold on
histogram(theta(:,end),[-5:5],'Normalization','pdf','FaceColor',colors.blue,'LineWidth',1.5);
[f, x] = ksdensity(theta(:,end));
mean_t = mean(theta(:,end));
sd_t = std(theta(:,end));
xline(mean_t,'LineWidth',4)
xline(mean_t-sd_t,'LineWidth',3,'LineStyle','--')
xline(mean_t+sd_t,'LineWidth',3,'LineStyle','--')
xline(mean_t-2*sd_t,'LineWidth',3,'LineStyle',':')
xline(mean_t+2*sd_t,'LineWidth',3,'LineStyle',':')
plot(x, f, 'LineWidth', 4, 'Color', colors.orange)
set(gca,'FontSize',22,'LineWidth',2)
[f, x] = ksdensity(theta(:,end));
%set(gca,'XTick',-5:5)
ylabel('PDE');
xlabel('SMART Score (θ)');
set(gca,'XTick',-5:5)
xlim([-5,5])
title('SMART Score','FontSize',26);
text(-0.19,1,'a','Units','normalized','FontSize',26,'FontWeight','bold')
text(0.02,0.98,...
    {'─   Mean', '╌ \pm1SD', '┈ \pm2SD'},...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top',...
    'FontSize',14);
box on 
grid on
hold off

nexttile;
hold on;
plot(1:length(mean_se),mean_se,'LineWidth', 5, 'Color', colors.blue)
plot(1:length(mean_rel),mean_rel,'LineWidth', 5, 'Color', colors.orange)
fill([1:length(mean_se), fliplr(1:length(mean_se))], [se_upper_bound, fliplr(se_lower_bound)], [0.2 0.6 0.8], ...
    'FaceAlpha', 0.5,'EdgeColor','none'); 
set(gca,'FontSize',22,'LineWidth',2)
fill([1:length(mean_rel), fliplr(1:length(mean_rel))], [rel_upper_bound, fliplr(rel_lower_bound)], [1 0.65 0.3], ...
    'FaceAlpha', 0.5,'EdgeColor','none'); 
xlim([2,20])
set(gca,'XTick',0:2:20)
ylabel('Estimate')
%ylim([-1.5,1.5])
%set(gca,'YTick',-1.5:0.25:1.5)
xlabel('Test Length')
title('SEM & Conditional Reliability','FontSize',26)
text(-0.15,1,'b','Units','normalized','FontSize',26,'FontWeight','bold')
%legend({'Standard Error', 'Reliability'},'FontSize',24)
lgd = legend('Standard Error', 'Reliability');
lgd.FontSize = 18;
drawnow;   % let MATLAB finalize legend size
lgd.Units = 'normalized';
ax = gca;
ax.Units = 'normalized';
axPos = ax.Position;
lgdPos = lgd.Position;
% Align SE corners
lgdPos(1) = axPos(1) + axPos(3) - lgdPos(3);
lgdPos(2) = axPos(2);
lgd.Position = lgdPos;
box on 
grid on
hold off
set(gcf, 'Units', 'pixels', 'Position', [1 1 1512 800])
exportgraphics(gcf, 'Figure_6.png', 'Resolution', 600)
%% Percentile ranks
ranking = [1 5 linspace(10,90,9) 95 99];
for i = 1:length(ranking)
    score_rank(i) = prctile(theta(:,end),ranking(i));
end
t_score = array2table([ranking;score_rank]','VariableNames',{'Percentile','Score'});
writetable(t_score,'percentile_scores.xlsx')
%% Correlation for different test length
idx = [];
for i=1:height(data)
    idx(i)=find(strcmpi(data.Row{i},userID));
end
theta = theta(idx,:);
vars_to_compare = {[6,15,16,18],[7:14],[2:5]};
significance = {[0.21,0.26,0.31],[0.189,0.232,0.279]};
figure
for k = 1:2
    rho_theta = NaN;
    for i =1:20
        n = 1;
        for v = vars_to_compare{k}
            rho_theta(i,n) = corr(theta(:,i),data{:,v},'rows','pairwise');
            n = n+1;
        end
    end
    subplot(1,2,k)
    plot(rho_theta,'LineWidth',5)
    ylabel('Correlation Coef.','FontSize',32);
    xlabel('Test Length','FontSize',24);
    set(gca,'FontSize',32,'LineWidth',2)
    xlim([1 size(rho_theta,1)])
    yline(significance{k}(1),'--',{'p = 0.05'})
    yline(significance{k}(2),'--',{'p = 0.01'})
    yline(significance{k}(3),'--',{'p = 0.001'})
    %yline(-0.21,'--',{'p = 0.05'})
    %yline(-0.26,'--',{'p = 0.01'})
    %yline(-0.31,'--',{'p = 0.001'})
    title('Correlations and Test Length','Fontsize',32,'Interpreter', 'none')
    legend(data.Properties.VariableNames(vars_to_compare{k}),'Location','northeastoutside',...
        'Fontsize',16,'Interpreter', 'none')
    box on
    grid on
end
%% Ability histogram
figure
subplot(1,2,1)
histogram(theta(:,end),[-5:5],'Normalization','pdf','FaceColor',colors.blue,'LineWidth',1.5);
set(gca,'FontSize',28,'LineWidth',2)
[f, x] = ksdensity(theta(:,end));
%set(gca,'XTick',-5:5)
ylabel('PDE');
xlabel('SMART Score (θ)');
mean_t = mean(theta(:,end));
%sd_l = mean_t - 2*std(theta(:,end));
%sd_h = mean_t + 2*std(theta(:,end));
xline(mean_t,'-',{'Mean Score'},'LineWidth',4,'FontSize',15)
%xline(sd_l,'-',{'2 SD'},'LineWidth',3)
%xline(sd_h,'-',{'2 SD'},'LineWidth',3)
title('Musical Affect Recognition Ability');
text(-0.2,1,'a','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
% Standard error
colors.blue = [0.1035, 0.4146, 0.6928];
colors.orange = [1 0.4980 0.0549,0.85];
for i = 1:20
    mean_dev(i) = nanmean(theta(:,end)-theta(:,i));
    std_dev(i) = nanstd(theta(:,end)-theta(:,i));
end
upper_bound = mean_dev + 2 * std_dev;
lower_bound = mean_dev - 2 * std_dev;
subplot(1,2,2)
hold on;
set(gca,'FontSize',28,'LineWidth',2)
fill([1:length(mean_dev), fliplr(1:length(mean_dev))], [upper_bound, fliplr(lower_bound)], colors.blue, ...
    'FaceAlpha', 0.5,'EdgeColor','none'); 
plot(1:length(mean_dev),mean_dev,'LineWidth', 5, 'Color', [0, 0, 0, 0.5])
xlim([2,19])
ylabel('Mean Deviation')
xlabel('Number of Trials')
title('Inter-Trial Ability Deviation')
text(-0.15,1,'b','Units','normalized','FontSize',30,'FontWeight','bold')
box on 
grid on
hold off
%% Percentage of data within -+ 1SD
upper_thx = nanmean(theta(:,end)) + nanstd(theta(:,end));
low_thx = nanmean(theta(:,end)) - nanstd(theta(:,end));
percentage = sum(theta(:,end)<upper_thx & theta(:,end)>low_thx)/size(theta,1);
%% Correct responses
figure
subplot(1,2,1)
histogram(mean(response,2),0:0.1:1,'Normalization','pdf','FaceColor',colors.blue,'LineWidth',1.5);
set(gca,'FontSize',28,'LineWidth',2)
[f, x] = ksdensity(mean(response,2));
%set(gca,'XTick',-5:5)
ylabel('PDE');
xlabel('Percentage');
mean_t = median(mean(response,2));
%sd_l = mean_t - 2*std(theta(:,end));
%sd_h = mean_t + 2*std(theta(:,end));
%xline(mean_t,'-',{'Mean Score'},'LineWidth',4,'FontSize',15)
%xline(sd_l,'-',{'2 SD'},'LineWidth',3)
%xline(sd_h,'-',{'2 SD'},'LineWidth',3)
title('Correct Responses');
text(-0.2,1,'a','Units','normalized','FontSize',30,'FontWeight','bold')
box on
grid on
