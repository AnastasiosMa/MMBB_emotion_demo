%% Run adaptive test simulation
trial_info = readtable('../data/output/binary_responses/final_trial_info.csv');
rasch_mirt = readtable('../data/output/binary_responses/irt_models/final_rasch_mirt.csv');
participant_scores = readtable('../data/output/binary_responses/irt_models/participantScores.csv');
emoNames = {'Angry','Fearful','Happy','Sad','Tender'};
item_difficulty = rasch_mirt{:,2};
item_emo = trial_info{:,2};
track = trial_info{:,end};
trialN = length(item_difficulty);
experiment = 1;
use_only_fixed_difficulty = 0;
starting_item_difficulty = [1 -1];
%% Create probability of correct and wrong sample answers
theta_step = 0.02;
theta_low = -6;
theta_high = 6;
theta_range = round(theta_low:theta_step:theta_high,2);
guessing = 0.5;
k=1;
for th = theta_range
    for j = 1:trialN
        p_correct(j,k) = guessing+(1-guessing)*(1/(1+exp(-(th-rasch_mirt{j,2}))));
        p_star(j,k) = 1/(1+exp(-(th-rasch_mirt{j,2})));
        p_incorrect(j,k) = 1-p_correct(j,k);
        product = p_correct(j,k)*p_incorrect(j,k);
        if k>1
            j_der1 = [p_correct(j,k)-p_correct(j,k-1)]./theta_step;
            information_test(j,k) = j_der1^2/product;
        end
    end
    k = k+1;
end
information_test(:,1) =  information_test(:,2);
%% Test Parameters
test_length = 30;
current_test_length = 20;
start_difficulty = 0;
optimizer = 1; %1 fixed difficulty, 2 ml optimizer
%simulation parameters
permutations = 100;
N = 1000;

%data matrices
if experiment==1
    th = 1.5*randn(permutations,N); %1000 participants, 1.5 std, mean 0
elseif experiment==2
    th = linspace(theta_low, theta_high,N);
    th = repmat(th,[permutations,1]);
end
track110 = nan(permutations,N,test_length);
th_est = nan(permutations,N,test_length);
trial_idx_selected = nan(permutations,N,test_length);
deviation_from_ground_truth = nan(permutations,N,test_length);
deviation_mistakes_only = nan(permutations,N,test_length);
emotion_history = nan(permutations,N,test_length);
optimizer_history = nan(permutations,N,test_length);
response_accuracy = nan(permutations,N);
deviation_from_item_difficulty = nan(permutations,N,test_length);
compare_information_vs_it = nan(permutations,N,test_length);
participant_responses = nan(permutations,N,test_length);
%% Test simulation
for p = 1:permutations
    disp(['Permutation: ' num2str(p)]);
    for k = 1: N
        %initialise parameters
        optimizer = 1;
        responses = nan(1,test_length);
        trial_idx = nan(1,test_length);
        epoch = 1;
        trial_idx = [];
        responses = [];
        track_selected = [];
        % crete emotion stratification
        emo_strat = [];
        for i = 1:test_length/5
            emo_strat = [emo_strat randperm(5)];
        end
        %find th_idx of TRUE ABILITY
        if th(p,k)<=theta_low
            th_idx = 1;
        elseif th(p,k)>=theta_high
            th_idx = length(theta_range);
        else
            th_idx = find(th(p,k)<=round(theta_range+theta_step,2) & th(p,k)>theta_range);
        end
        %set th_est =0
        th_est_idx = round(length(theta_range)/2);
        th_est(p,k,epoch) = 0;
        starting_item_difficulty = starting_item_difficulty(randperm(length(starting_item_difficulty)));
        while epoch<=test_length
            if epoch>length(starting_item_difficulty)
                %remove trials of previously listened excerpts
                trial_candidates = find(item_emo==emo_strat(epoch) & all(track_selected~=track,2));
                %shuffle candidate trials
                trial_candidates = trial_candidates(randperm(length(trial_candidates)));

                %[~,min_dist] = sort(abs(th_est(p,k,epoch-1) - item_difficulty(trial_candidates)));
                %trial_idx(epoch) = trial_candidates(min_dist(1));
                %temp_trial = trial_idx(epoch);
                [~,trial_idx(epoch)] = max(information_test(trial_candidates,th_est_idx));
                trial_idx(epoch) = trial_candidates(trial_idx(epoch));
                %compare_information_vs_it(p,k,epoch) = item_difficulty(temp_trial)-item_difficulty(trial_idx(epoch));
            else
                if epoch==1
                    trial_candidates = find(item_emo==emo_strat(epoch));
                else
                    trial_candidates = find(item_emo==emo_strat(epoch) & all(track_selected~=track,2));
                end
                trial_candidates = trial_candidates(randperm(length(trial_candidates)));
                [~,min_dist] = sort(abs(starting_item_difficulty(epoch) - item_difficulty(trial_candidates)));
                trial_idx(epoch) = trial_candidates(min_dist(randi(5,1)));
            end
            deviation_from_item_difficulty(p,k,epoch) = th(p,k)-item_difficulty(trial_idx(epoch));
            track_selected(epoch) = track(trial_idx(epoch));
            %propabilities of participant response
            responses(epoch) = randsrc(1,1,[0,1;p_incorrect(trial_idx(epoch),th_idx),...
                p_correct(trial_idx(epoch),th_idx)]);
            if responses(epoch)
                deviation_from_ground_truth(p,k,epoch) = p_incorrect(trial_idx(epoch),th_idx);
            else
                deviation_from_ground_truth(p,k,epoch) = -p_correct(trial_idx(epoch),th_idx);
            end
            %check for optimizer change
            if all(responses) | all(~responses)
                optimizer = 1;
            elseif ~use_only_fixed_difficulty
                optimizer = 2;
            end
            if epoch>1
                if isempty(th_est(p,k,epoch))
                    keyboard
                end
                [th_est(p,k,epoch), th_est_idx] = ml_optimizer(th_est(p,k,epoch-1),optimizer,...
                    responses,trial_idx,p_star,p_incorrect,th(p,k));
            else
                if isempty(th_est(p,k,epoch))
                    keyboard
                end
                [th_est(p,k,epoch), th_est_idx] = ml_optimizer(th_est(p,k,epoch),optimizer,...
                    responses,trial_idx,p_star,p_incorrect,th(p,k));
            end
            if isnan(th_est(p,k,epoch)) | isempty(th_est_idx)
                keyboard
            end
            optimizer_history(p,k,epoch) = optimizer;
            epoch = epoch+1;
        end
        response_accuracy(p,k) = mean(responses);
        participant_responses(p,k,:) = responses;
        emotion_history(p,k,:) = emo_strat;
        trial_idx_selected(p,k,:) = trial_idx;
        track110(p,k,:) = track_selected;
    end
end
%% Compare optimizer with Aggregate Difficulty Scores
difficulty_scores = rescale(item_difficulty);

for p = 1:permutations
    for k = 1:N
        for e = 1:test_length
            correct_score = 0;
            incorrect_score = 0;
            if ~isempty(find(participant_responses(p,k,1:e)))
                correct_idx = find(participant_responses(p,k,1:e));
                correct_score = sum(difficulty_scores(trial_idx_selected(p,k,correct_idx)));
            end
            if ~isempty(find(participant_responses(p,k,1:e)==0))
                incorrect_idx = find(participant_responses(p,k,1:e)==0);
                incorrect_score = sum(difficulty_scores(trial_idx_selected(p,k,incorrect_idx)));
                ads(p,k,e) = (correct_score-incorrect_score)/e;
            end
        end
    end
end

%get mean error and corr of th and th_est
for i = 1:test_length
    th_model = th_est(:,:,i);
    ads_model = ads(:,:,i);
    rho(i) = corr(th_model(:),th(:),'rows','pairwise');
    rho_ads(i) = corr(ads_model(:),th(:),'rows','pairwise');
    mae(:,i) = th(:) - th_model(:);
    deviation = mean(deviation_from_ground_truth(:,:,1:i),3);
    rho_div(i) = corr(deviation(:),th_model(:),'rows','pairwise');
    dev_mean(i) = mean(mean(deviation));
    [b(i,:),~,~,~,stats] = regress(rescale(th(:)),[rescale(th_model(:)),rescale(deviation(:)),ones(numel(th_model),1)]);
    %[b(i,:),~,~,~,stats] = regress(rescale(th(:)),[rescale(th_model(:)),deviation_low(:),ones(numel(th_model),1)]);
    r_sq_dev(i) = stats(1);
    item_dev(i) = mean(mean(deviation_from_item_difficulty(:,:,i)));
end
%% Plots
%distribution of generated abilities
figure
subplot(1,2,1)
histogram(th(:),theta_low:theta_high)
set(gca,'FontSize',32,'LineWidth',2)
title('Generated participants ability')
box on
xlabel('Θ')
grid on
hold off

subplot(1,2,2)
histogram(th_est(:,:,current_test_length),theta_low:theta_high)
set(gca,'FontSize',32,'LineWidth',2)
title('Optimizer participants ability')
box on
xlabel('Θ')
grid on
hold off

figure
deviation = mean(deviation_from_ground_truth(:,:,1:current_test_length),3);
histogram(deviation(:))
set(gca,'FontSize',32,'LineWidth',2)
title('Deviation from ground truth')
box on
xlabel('Deviation')
ylabel('Count')
grid on
hold off

figure
subplot(1,2,1)
hold on
plot(rho,'LineWidth',5)
ylabel('Correlation Coefficient','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
xline(current_test_length,'LineWidth',5)
title('True ability and ML Optimizer')
box on
grid on
hold off

subplot(1,2,2)
hold on
plot(rho_ads,'LineWidth',5)
ylabel('Correlation Coefficient','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
xline(current_test_length,'LineWidth',5)
title('True ability and Aggregate Difficulty Scores')
box on
grid on
hold off

figure
hold on
plot(r_sq_dev,'LineWidth',5)
plot(b(:,1),'LineWidth',5)
ylabel('Estimate','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
xline(current_test_length,'LineWidth',5)
title('R Square: θ(sim) = b1*θ(opt) + b2*deviation')
legend('R Square','b1','Location','best')
box on
grid on
hold off

figure
histogram(mae(:,current_test_length),theta_low:theta_high)
set(gca,'FontSize',32,'LineWidth',2)
title(['Mean Error of θ-θhat for Trial Length = ', num2str(current_test_length)])
box on
xlabel('Error','FontSize',32)
grid on
hold off

figure
hold on
plot(mean(mae),'LineWidth',5)
ylabel('Mean error','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
xline(current_test_length,'LineWidth',5)
title('Mean error (θ-θ hat)')
box on
grid on
hold off

figure
hold on
scatter(th_est(1,:,current_test_length),th(1,:),100,'filled');
xlabel('Optimizer estimates','FontSize',24);
ylabel('True estimates','FontSize',24);
title({['Scatterplot of 1 permutation, TestLength=20'], ['r=' num2str(round(corr(th_est(1,:,current_test_length)',...
    th(1,:)','type','Spearman'),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
h = lsline();
h.LineWidth = 5;h.Color = 'k';
%xlim([-5 5]);
%ylim([-5 5])
box on
grid on
hold off

%optimizer ratio
for i =1:test_length
    op_data = optimizer_history(:,:,i);
    op_mean(i)=mean(op_data(:)-1);
end

figure
hold on
plot(op_mean,'LineWidth',5)
ylabel('Optimizer ratio','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 length(op_mean)])
xline(current_test_length,'LineWidth',5)
title('Optimizer/Fixed Difficulty')
box on
grid on
hold off

figure
hold on
plot(item_dev,'LineWidth',5)
ylabel('Deviation','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
xline(current_test_length,'LineWidth',5)
title('Mean θ-b across Test Length')
box on
grid on
hold off

if experiment==2
    for i =15:test_length
        item_dif_dev{i} = mean(deviation_from_item_difficulty(:,:,15:i),3);
    end
    figure
    plot(mean(item_dif_dev{current_test_length}),'LineWidth',5)
    ylabel('Mean Deviation','FontSize',32);
    xlabel('Θ','FontSize',24);
    set(gca,'FontSize',32,'LineWidth',2)
    title('Item Difficulty Deviation ')
    box on
    grid on
    hold off

    %plot item difficulty selected
    diff_selected = item_difficulty(trial_idx_selected);

    figure
    plot(mean(mean(diff_selected(:,:,10:current_test_length),3)),'LineWidth',5)
    ylabel('Mean Item Difficulty','FontSize',32);
    xlabel('Θ','FontSize',24);
    set(gca,'FontSize',32,'LineWidth',2)
    title('Item Difficulty across Theta')
    box on
    grid on
    hold off
end
%% Check results for different levels of theta
th_model = th_est(:,:,current_test_length);
th_model = th_model(:);
th_temp = th(:);
i = 1;
for k = theta_low:theta_high
    idx = find(th_temp>k & th_temp<k+1);
    sd_theta(i) = std(th_temp(idx) - th_model(idx));
    mean_theta(i) = mean(th_temp(idx) - th_model(idx));
    groupN(i) = length(idx);
    i = i+1;
end

figure
plot([mean_theta',sd_theta'],'LineWidth',5)
ylabel('Metric','FontSize',32);
xlabel('Theta','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2,'XTick',1:length(theta_low:theta_high),...
    'XTickLabel',theta_low:theta_high)
title('Mean and SD for different theta values')
box on
grid on
legend({'Mean','SD'},'Location','best')

%deviation across theta
deviation = deviation_from_ground_truth(:,:,current_test_length);
deviation = deviation(:);
th_temp = th(:);
i = 1;
for k = theta_low:theta_high
    idx = find(th_temp>k & th_temp<k+1);
    positive_deviation(i) = sum(deviation(idx)>0)/length(idx);
    negative_deviation(i) = sum(deviation(idx)<-0.5)/length(idx);
    groupN(i) = length(idx);
    i = i+1;
end

figure
plot([positive_deviation',negative_deviation'],'LineWidth',5)
ylabel('Metric','FontSize',32);
xlabel('Deviation from ground truth','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2,'XTick',1:length(theta_low:theta_high),...
    'XTickLabel',theta_low:theta_high)
title('Positive and negative deviationfor different theta values')
box on
grid on
legend({'Positive Dev','Negative Dev'},'Location','best')
%% Check big mistakes
p = 1; %choose a permutation
th_model = th_est(p,:,current_test_length);
[~,i] = sort(th(p,:)-th_model(:)'); %find biggest mistakes;
th(p,i(1))-th_model(p,i(1))
%data from biggest mistake
p_resp = participant_responses(p,i(1),1:current_test_length);
trial_id = trial_idx_selected(p,i(1),1:current_test_length);
t = th_est(p,i(1),1:current_test_length);
if th(p,i(1))<=theta_low
    th_idx = 1;
elseif th(p,i(1))>=theta_high
    th_idx = length(theta_range);
else
    th_idx = find(th(p,i(1))<=theta_range+theta_step & th(p,i(1))>theta_range);
end
d = [item_difficulty(trial_id(:)),trial_id(:),p_resp(:),t(:),p_correct(trial_id(:),th_idx)];
th(p,i(1))
%check deviations
tl = 2; %test length
p_resp = squeeze(participant_responses(p,:,1:tl));
i = find(p_resp(:,1) == 1 & p_resp(:,2) == 0);
trial_id = trial_idx_selected(p,i(1),1:tl);
dev = deviation_from_ground_truth(p,i(1),1:tl);
t = th_est(p,i(1),1:tl);
if th(p,i(1))<=theta_low
    th_idx = 1;
elseif th(p,i(1))>=theta_high
    th_idx = length(theta_range);
else
    th_idx = find(th(p,i(1))<=theta_range+theta_step & th(p,i(1))>theta_range);
end
d = [dev(:),item_difficulty(trial_id(:)),trial_id(:),p_resp(i(1),1:tl)',t(:),p_correct(trial_id(:),th_idx)];
th(p,i(1))
%% Check distribution (number of appearances) of all items
items_selected = trial_idx_selected(:,:,1:current_test_length);
items_selected = items_selected(:);
[item_frequency, gr] = groupcounts(items_selected);

m = max(groupcounts(item_emo));
gc = groupcounts(item_emo)/m;

figure
hold on
scatter(item_difficulty(gr),item_frequency,100,'filled');
set(gca,'FontSize',32,'LineWidth',2)
xlabel('Item Difficulty');
ylabel('Number of occurences');
title('Occurences per item')
yline(mean(item_frequency),'-',{'Mean Occurence'},'LineWidth',5)
yline(median(item_frequency),'-',{'Median Occurence'},'LineWidth',5)
box on
grid on

figure
gscatter(item_difficulty(gr),item_frequency,item_emo(gr),lines(length(unique(item_emo))),'..',35);
set(gca,'FontSize',32,'LineWidth',2)
xlabel('Item Difficulty');
ylabel('Number of occurences');
title('Occurences per item')
yline(mean(item_frequency),'-',{'Mean Occurence'},'LineWidth',5)
yline(median(item_frequency),'-',{'Median Occurence'},'LineWidth',5)
legend(emoNames,'Location','best')
box on
grid on

figure
gscatter(item_difficulty(gr),item_frequency.*gc(item_emo(gr)),item_emo(gr),lines(length(unique(item_emo))),'..',35);
set(gca,'FontSize',32,'LineWidth',2)
xlabel('Item Difficulty');
ylabel('Number of occurences');
title('Occurences per item (Adjusted per emotion)')
yline(mean(item_frequency),'-',{'Mean Occurence'},'LineWidth',5)
yline(median(item_frequency),'-',{'Median Occurence'},'LineWidth',5)
legend(emoNames,'Location','best')
box on
grid on