% HARDCODED MICROPROMS-CHECK DATA FILE
%publish('validation_analysis.m','format','pdf','showCode',false);
warning('off')
qualtrics_data = readtable('data/qualtrics_results.csv');
qualtrics2 = readtable('qualtrics/preprocessed_ERI_part2_combined.csv');
microproms_data = readtable('data/MICROPROMS_results.csv');
prolific_demographics = readtable('data/prolific_export.csv');
jspsych_dir = 'data/jatos_results/';
jatos_metadata = 'data/jspsych_metadata.json';
pilot_participants = [23:31,33,34];
duplicate_participants = {'60fe9a0daa398bbc0d610dc0'};
pilot_ids = {'662fc8f6aa3d3b206672e28d','5c181b8d905b250001be42ca','66a27a561b8ed92d82f3b97a',...
    '5ea1a2e7939362055f0c325a','5dc9c2c0b9fad36ddca37632','663a3ee80ab62f61353a10f0','562a00dac8ffc20012513fbe',...
    '55d22025cc2b18000c0b9d9c','6715521230e88bbd98ef83ba','62d6d87d91cf24a3c1c0f10c','5d49d17b3dad1f0001e2aba1'};
tas_describe_feelings_idx = [2, 4, 11, 12, 17];
tas_identify_feelings_idx = [1, 3, 6, 7, 9, 13, 14];
tas_external_thinking_idx = [5, 8, 10, 15, 16, 18, 19, 20];
tas_negative_questions = [4, 5, 10, 18, 19];
msi_negative_questions = [1,2];
% Filter qualtrics responses
for i = 1:height(qualtrics_data)
    q_idx(i)=length(qualtrics_data{i,'UserID'}{1})>10;
end
qualtrics_data = qualtrics_data(q_idx,:);
for i = 1:length(duplicate_participants)
if any(strcmpi(qualtrics_data.UserID,duplicate_participants{i}))
   qualtrics_data.UserID{find(strcmpi(qualtrics_data.UserID,duplicate_participants{i}),1)}=[duplicate_participants{i},'_dble']; 
end
end
qualtrics_data.Row = qualtrics_data.UserID;
qualtrics_data = qualtrics_data(:,{'score_facial','score_vocal'});
qualtrics2.Row = qualtrics2.UseriD;
qualtrics2 = qualtrics2(:,{'score_facial','score_vocal'});
qualtrics_data = vertcat(qualtrics_data,qualtrics2);
%correlate face_vocal factors
rho = corr(qualtrics_data{:,'score_facial'},qualtrics_data{:,'score_vocal'},'rows','pairwise');
%disp(['Correlation between facial and vocal: ',num2str(round(rho,2))])
% Filter MICROPROMS responses
for i = 1:height(microproms_data)
    m_idx(i)=length(microproms_data{i,1}{1})>10;
end
microproms_data = microproms_data(m_idx,:);
for i = 1:length(duplicate_participants)
if any(strcmpi(microproms_data.Var1,duplicate_participants{i}))
   microproms_data.Var1{find(strcmpi(microproms_data.Var1,duplicate_participants{i}),1)}=[duplicate_participants{i},'_dble']; 
end
end
microproms_data.Row = microproms_data.Var1;
%microproms_data.x_time___CaliStart_ = [];
microproms_data.Var1 = [];
microproms_col_names = {'mProms_total','mProms_melody','mProms_timbre','mProms_tempo',...
                        'mProms_beat','mProms_tuning','mProms_rhythm','mProms_pitch'}; %HARDCODED
microproms_data.Properties.VariableNames = microproms_col_names;
% Demographic data
%prolific_demographics = prolific_demographics(strcmpi('APPROVED',prolific_demographics.Status),:);
for i = 1:height(prolific_demographics)
    if strcmpi(prolific_demographics{i,'Sex'},'Male')
       prolific_demographics{i,'Sex_bin'} = 1;
    elseif strcmpi(prolific_demographics{i,'Sex'},'Female')
       prolific_demographics{i,'Sex_bin'} = 0;
    else
       prolific_demographics{i,'Sex_bin'} = 2;
    end
end
prolific_demographics.Row = prolific_demographics.ParticipantId;
prolific_vars = {'TimeTaken','TotalApprovals','Age','Sex_bin'};
prolific_demographics = prolific_demographics(:,prolific_vars);
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
        if t ==trial_idx(2)
            msi(i,1) = jspsych_data.trials{t}.response.MSI1;
            msi(i,2) = jspsych_data.trials{t}.response.MSI2;
        elseif t ==trial_idx(3)
            msi_names = fieldnames(jspsych_data.trials{t}.response);
            for n = 1:length(msi_names)
                if ~isempty(str2num(getfield(jspsych_data.trials{t}.response,msi_names{n})))
                    msi(i,n+2) = str2num(getfield(jspsych_data.trials{t}.response,msi_names{n}));
                else
                    msi(i,n+2) = nan;
                end
            end
        elseif t ==trial_idx(4)
            tas_names = fieldnames(jspsych_data.trials{t}.response);
            for n = 1:length(tas_names)
                if ~isempty(getfield(jspsych_data.trials{t}.response,tas_names{n}))
                    tas(i,n) = getfield(jspsych_data.trials{t}.response,tas_names{n});
                else
                    tas(i,n) = nan;
                end
            end
        elseif t ==trial_idx(5)
            tas_names = fieldnames(jspsych_data.trials{t}.response);
            for n = 1:length(tas_names)
                if ~isempty(getfield(jspsych_data.trials{t}.response,tas_names{n}))
                    tas(i,n+10) = getfield(jspsych_data.trials{t}.response,tas_names{n});
                else
                    tas(i,n+10) = nan;
                end
            end
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
disp(['Test time (minutes), Median: ',num2str(round(nanmedian(prolific_demographics.TimeTaken/60),2))])
%disp(['Test time, MAD: ',num2str(round(mad(prolific_demographics.TimeTaken/60),2))])
%disp(['Exclusion Threshold(3 MADs): ',num2str(round(nanmedian(prolific_demographics.TimeTaken/60)-mad(prolific_demographics.TimeTaken/60)*3,2))])

%create matrix of thetas
for k = 1:(length(jspsych_table))
    theta(k,:) = jspsych_table{k}.Theta;
    difficulty(k,:) = jspsych_table{k}.Difficulty;
    response(k,:) = jspsych_table{k}.response;
end

%Drop non prolific cases
theta = theta(valid_id,:);
msi = msi(valid_id,:);
tas = tas(valid_id,:);
userID = userID(valid_id);

%tas convert negatively keyed terms
tas(:,tas_negative_questions) = 4 - tas(:,tas_negative_questions);
tas = tas+1;
%tas create factor scores
tas_factors(:,1) = nansum(tas(:,tas_describe_feelings_idx),2);
tas_factors(:,2) = nansum(tas(:,tas_identify_feelings_idx),2);
tas_factors(:,3) = nansum(tas(:,tas_external_thinking_idx),2);
tas_factors(:,4) = nansum(tas,2);

%convert starttime
date = datetime(starttime/1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
date = date(valid_id);
week_day = weekday(date);
hour_num = hour(date);
for i = 1:length(hour_num)
    if hour_num(i)<12
       timeofday{i} = 'Morning';
    elseif hour_num(i)>=12 & hour_num(i)<15
       timeofday{i} = 'Noon';
    elseif hour_num(i)>=15 & hour_num(i)<18
       timeofday{i} = 'Afternoon';
    elseif hour_num(i)>=18
       timeofday{i} = 'Evening';
    end
end
% Construct msi factors
msi_scoring = {[0;1;2;3;4;6;10],[0;1;2;3;4;5;6],[0;0.5;1;1.5;2;3;5],[0;0.5;1;2;3;4;7],...
            [0;0.5;1;2;3;6;10]};
for k = 1:size(msi,1)
for i = 1:5
    msi(k,i+2) = sum(msi(k,i+2)>=msi_scoring{i})-1;
    if msi(k,i+2) ==-1
       msi(k,i+2) = NaN; 
    end
end
end
msi(:,msi_negative_questions) = 6 - msi(:,msi_negative_questions);
msi = msi+1;
msi_factor = nansum(msi,2);
disp('Ability level (Theta) Summary statistics')
disp(['Median Theta: ', num2str(round(median(theta(:,end)),2))])
disp(['Std Theta: ', num2str(round(std(theta(:,end)),2))])
disp(['Max Theta: ', num2str(round(max(theta(:,end)),2))])
disp(['Min Theta: ', num2str(round(min(theta(:,end)),2))])

% Concatenate factors
%jspsych
jspsych = array2table([theta(:,end),tas_factors,msi_factor],'VariableNames',{'Theta',...
    'Tas_DescEm','Tas_IdEm','Tas_ExtTh','Tas_Total','MSI_MTraining'},'RowNames',userID);
%microproms
microproms_sorted = array2table(nan(size(jspsych,1),size(microproms_data,2)),'VariableNames',microproms_data.Properties.VariableNames);
qualtrics_sorted = array2table(nan(size(jspsych,1),size(qualtrics_data,2)),'VariableNames',qualtrics_data.Properties.VariableNames);
demographics_sorted = array2table(nan(size(jspsych,1),size(prolific_demographics,2)),'VariableNames',prolific_demographics.Properties.VariableNames);
for i = 1:length(jspsych.Row)
    if find(strcmpi(jspsych.Row{i},microproms_data.Row))
       microproms_sorted(i,:) = microproms_data(find(strcmpi(jspsych.Row{i},microproms_data.Row)),:); 
    end
    if find(strcmpi(jspsych.Row{i},qualtrics_data.Row))
       qualtrics_sorted(i,:) = qualtrics_data(find(strcmpi(jspsych.Row{i},qualtrics_data.Row)),:); 
    end
    if find(strcmpi(jspsych.Row{i},prolific_demographics.Row))
       demographics_sorted(i,:) = prolific_demographics(find(strcmpi(jspsych.Row{i},prolific_demographics.Row)),:); 
    end
end
data = [jspsych,microproms_sorted,qualtrics_sorted,demographics_sorted];
%remove missing values
data.TotalApprovals = [];
%% Demographics
disp(['N: ', num2str(height(jspsych))])
disp(['N (complete data): ', num2str(sum(~isnan(qualtrics_data.score_vocal)))])
disp(['Mean Age: ', num2str(round(nanmean(prolific_demographics.Age),2)), ...
    '  Std: ', num2str(round(nanstd(prolific_demographics.Age),2))])
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
title('SEM & Conditional Reliability','FontSize',30)
text(-0.15,1,'b','Units','normalized','FontSize',30,'FontWeight','bold')
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
%% Plot test time
predictors = [5,6,7,15,16,18];
cats = {'Morning','Noon','Afternoon','Evening'};
for i = 1:4
    time_means(i,:)=nanmean(table2array(data(find(strcmpi(timeofday,cats{i})),[1,predictors])));
end
time_means(:,1) = time_means(:,1)*10; %multiply theta by 10 for plotting
figure
plot(time_means,'LineWidth',5)
ylabel('Mean Scores','FontSize',32);
xlabel('Time of Day','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
set(gca,'XTick',1:4,'XTickLabel',cats,'FontSize',18)
title('Time & Test Scores','Fontsize',32,'Interpreter', 'none')
legend(data.Properties.VariableNames([1,predictors]),'Location','northeastoutside',...
    'Fontsize',16,'Interpreter', 'none')
box on
grid on

%day of the week
for i = 1:5
    day_means(i,:)=nanmean(table2array(data(find(week_day==i+1),[1,predictors])));
end
day_means(:,1) = day_means(:,1)*10; %multiply theta by 10 for plotting
figure
plot(day_means,'LineWidth',5)
ylabel('Mean Scores','FontSize',32);
xlabel('Day of the week','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
set(gca,'XTick',1:5,'XTickLabel',{'Mon','Tue','Wed','Thu','Fri'},'FontSize',18)
title('Day & Test Scores','Fontsize',32,'Interpreter', 'none')
legend(data.Properties.VariableNames([1,predictors]),'Location','northeastoutside',...
    'Fontsize',16,'Interpreter', 'none')
box on
grid on
%% Make correlation table
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
writetable(corr_tbl,'../paper/corr_table.csv')
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
