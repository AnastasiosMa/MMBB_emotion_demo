data = readtable('../data/output/binary_responses/final_binary_responses.csv');
covariates = data(:,end-2:end);
data = data(:,1:end-3);
trial_info = readtable('../data/output/binary_responses/final_trial_info.csv');
trialN = size(data,2);
emoNames = {'Angry','Fearful','Happy','Sad','Tender'};
difficulty_item = mean(data{:,1:trialN},'omitnan');
rasch_mirt = readtable('../data/output/binary_responses/irt_models/final_rasch_mirt.csv');
participant_scores = readtable('../data/output/binary_responses/irt_models/participantScores.csv');
guessing = 0.5;

%HARDCODED LOCATION OF PARTICIPANTS
fi_groups = [repmat(1,1,67),repmat(2,1,47)];
spa_groups = repmat(3,1,length(115:243));
spa_groups(find(isnan(data{115:243,1}))) = 4;
groups = [fi_groups,spa_groups];

[~,ia]  = unique(groups);
group_binary = ~isnan(data{ia,1:trialN});
max_item_N = max(sum(group_binary,2));
%% Optimizer denominator
%create probability of correct and wrong sample answers
theta_step = 0.05;
theta_low = -8;
theta_high = 8;
theta_range = round(theta_low:theta_step:theta_high,2);
k=1;
for th = theta_range
    for j = 1:trialN
        p_correct(j,k) = guessing+(1-guessing)*(1/(1+exp(-(th-rasch_mirt{j,2}))));
        p_incorrect(j,k) = 1-p_correct(j,k);
        %information_test(j,k) = p_correct(j,k)*p_incorrect(j,k);
        %if k>2
        %    j_der1 = [p_correct(j,k)-p_correct(j,k-1),p_correct(j,k-1)-p_correct(j,k-2)]./...
        %        [theta_step,theta_step];
        %    j_der2 = (j_der1(1)-j_der1(2))/theta_step;
        %    j_final(j,k) = (j_der1(1)*j_der2)/information_test(j,k);
        %end
    end
    k = k+1;
end
%j_final(:,1) =  j_final(:,3); j_final(:,2) =  j_final(:,3);
%th = th_idx, j = item order, u = participant responses th = theta value
%index optimizer formula
optimizer = @(th,u,j) sum(u-p_correct(j,th)')/(-sum(p_correct(j,th)'.*p_incorrect(j,th)'));
%bias_correction = @(th,j) sum(j_final(j,th))/(2*sum(information_test(j,th)));
%% Test optimizer
init_th = 0; %initial starting point of theta
learning_rate = 0.5;
min_criterion = theta_step; %optimizer stops if Δθ < criterion
n_iter = 500; %number of iterations before manual stop
%timings = nan(size(data,1),size(data,2));
iterations = nan(max_item_N,n_iter);
optimizer_training_history = cell(size(data,1),500);
optimizer_weight_adjustment_history = cell(size(data,1),500);
for permutations = 1:500
    disp(permutations)
    perms = randperm(size(data,2));
    for participant = 1:size(data,1) %per participant
        optimizer_training_history{participant,permutations} = iterations;
        participant_data = data{participant,perms};
        p_resp  = find(~isnan(participant_data));%non-empty responses for ther participant
        for trial = 2:length(p_resp) %per participant trial
            delta_th = 1;
            th = init_th;
            iter = 1;
            tStart = tic;
            while abs(delta_th)>min_criterion
                if iter>n_iter
                    break
                end
                if th<theta_low
                    th_idx = 1;
                    th = theta_low;
                    if ~any(participant_data(p_resp(1:trial))) %if all answers wrong, break
                        break
                    end
                elseif th>theta_high
                    th_idx = length(theta_range);
                    th = theta_high;
                    if all(participant_data(p_resp(1:trial)))%if all answers correct, break
                        break
                    end
                else
                    %find discrete index for theta within range
                    th_idx = find(th<=theta_range+theta_step & th>theta_range);
                end
                if isempty(th_idx)
                    keyboard
                end
                delta_th = optimizer(th_idx,participant_data(p_resp(1:trial)),...
                    p_resp(1:trial));
                th = th-delta_th;
                optimizer_training_history{participant,permutations}(trial-1,iter) = th;
                iter = iter+1;
            end
            %timings(participant,trial-1) = toc(tStart);
            %optimizer_weight_adjustment_history{participant,permutations}(trial-1) = th + bias_correction(th_idx,p_resp(1:trial));
        end
    end
end
%% Plots
%get theta estimates for all responses of each participant
for k = 1:size(optimizer_training_history,1)
    row = max(find(~isnan(optimizer_training_history{k,1}(:,1))));
    col = max(find(~isnan(optimizer_training_history{k,1}(row,:))));
    theta_end(k) = optimizer_training_history{k,1}(row,col);
    %weighted_theta_end(k) = optimizer_weight_adjustment_history{k,1}(end);
end
difficulty_p = mean(data{:,1:trialN}','omitnan');
figure
scatter(difficulty_p,theta_end,100,'filled')
rho = corr(difficulty_p',theta_end');
xlabel('Participant ability (response accuracy)',...
    'FontSize',24)
ylabel('Participant ability (optimizer estimate)','FontSize',24)
title(['Optimizer and ground truth correlation: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
box on
grid on

%nationality
figure
gscatter(difficulty_p,theta_end,covariates{:,1},lines(2),'..',35)
rho = corr(difficulty_p',theta_end');
xlabel('Participant ability (response accuracy)',...
    'FontSize',24)
ylabel('Participant ability (optimizer estimate)','FontSize',24)
title(['Optimizer and ground truth correlation: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
legend({'Finnish','Spanish'},'Location','best')
box on
grid on

%groups
figure
gscatter(difficulty_p,theta_end,groups,lines(4),'..',35)
rho = corr(difficulty_p',theta_end');
xlabel('Participant ability (response accuracy)',...
    'FontSize',24)
ylabel('Participant ability (optimizer estimate)','FontSize',24)
title(['Participant groups: r = ', num2str(round(rho,2))],...
    'FontSize',30)
set(gca,'FontSize',32,'LineWidth',2)
legend({'Finns Group 2','Finns Group 2','Spanish Group 1','Spanish Group 2'},'Location','best')
box on
grid on
%% Plot performance across test Length
%correlations with random permutations
for p = 1:500
    disp(p)
    for i = 1:max_item_N
        for k = 1:size(optimizer_training_history,1)
            col = max(find(~isnan(optimizer_training_history{k,p}(i,:))));
            if ~isempty(optimizer_training_history{k,p}(i,col))
                theta(k,i) = optimizer_training_history{k,p}(i,col);
                %weighted_theta(k,i) = optimizer_weight_adjustment_history{k,p}(i);
            else
                theta(k,i) = NaN;
                %weighted_theta(k,i) = NaN;
            end
        end
        rho(p,i) = corr(difficulty_p',theta(:,i),'rows','pairwise');
        %rho_weighted(p,i) = corr(difficulty_p',weighted_theta(:,i),'rows','pairwise');
    end
end

figure
hold on
plot(mean(rho)','LineWidth',5)
ylabel('Correlation Coefficient','FontSize',32);
xlabel('Number of Items in Theta estimation','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
title('Mean r across permutations')
low_ci = mean(rho)-(1.5*std(rho));
high_ci = mean(rho)+(1.5*std(rho));
%fill([1:size(rho,2) fliplr(1:size(rho,2))], ...
 %   [low_ci fliplr(high_ci)], 'red', 'FaceColor','r','FaceAlpha',0.3,'LineWidth',1)
box on
grid on
hold off
%% Correlation of optimizer scores and model ability estimates
figure
hold on
scatter(theta_end',participant_scores{:,1},100,'filled');
xlabel('Optimizer estimates','FontSize',24);
ylabel('Model estimates','FontSize',24);
title({['Optimizer and Model ability estimates'], ['r=' num2str(round(corr(theta_end',...
    participant_scores{:,1}),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
h = lsline();
h.LineWidth = 5;h.Color = 'k';
xlim([-3 3]);
ylim([-3 3])
box on
grid on
hold off