%% Run adaptive test simulation
trial_info = readtable('../data/output/binary_responses/final_trial_info.csv');
rasch_mirt = readtable('../data/output/binary_responses/irt_models/final_rasch_mirt.csv');
item_difficulty = rasch_mirt{:,2};
item_emo = trial_info{:,2};
track = trial_info{:,end};
trialN = length(item_difficulty);
%% Create probability of correct and wrong sample answers
theta_step = 0.05;
theta_low = -6;
theta_high = 6;
theta_range = round(theta_low:theta_step:theta_high,2);
guessing = 0.5;
k=1;
for th = theta_range
    for j = 1:trialN
        p_correct(j,k) = guessing+(1-guessing)*(1/(1+exp(-(th-rasch_mirt{j,2}))));
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
start_difficulty = 0;
optimizer = 1; %1 fixed difficulty, 2 ml optimizer
%simulation parameters
permutations = 1000;
N = 1000;

%data matrices
th = 1.5.*randn(N,permutations); %1000 participants, 1.5 std, mean 0
track110 = nan(permutations,N,test_length);
th_est = nan(permutations,N,test_length);
trial_idx_selected = nan(permutations,N,test_length);
emotion_history = nan(permutations,N,test_length);
optimizer_history = nan(permutations,N,test_length);
response_accuracy = nan(permutations,N);
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
        %find start trial at difficulty 0
        [~,start_trial_candidates] = sort(abs(item_difficulty(find(item_emo==emo_strat(1)))));
        trial_idx(epoch) = start_trial_candidates(randi(4,1));
        %find th_idx of TRUE ABILITY
        if th(p,k)<theta_low
            th_idx = 1;
        elseif th(p,k)>theta_high
            th_idx = length(theta_range);
        else
            th_idx = find(th(p,k)<=theta_range+theta_step & th(p,k)>theta_range);
        end
        while epoch<=test_length
            if epoch>1
                %remove trials of previously listened excerpts
                trial_candidates = find(item_emo==emo_strat(epoch) & all(track_selected~=track,2));
                %shuffle candidate trials
                trial_candidates = trial_candidates(randperm(length(trial_candidates)));
                [~,trial_idx(epoch)] = max(information_test(trial_candidates,th_est_idx));
                trial_idx(epoch) = trial_candidates(trial_idx(epoch));
            else
                th_est_idx = round(length(theta_range)/2);
                th_est(p,k,epoch) = 0;
            end
            track_selected(epoch) = track(trial_idx(epoch));
            %propabilities of participant response
            responses(epoch) = randsrc(1,1,[0,1;p_incorrect(trial_idx(epoch),th_idx),...
                p_correct(trial_idx(epoch),th_idx)]);
            %check for optimizer change
            if all(responses) | all(~responses)
                optimizer = 1;
            else
                optimizer = 2;
            end
            if epoch>1
                [th_est(p,k,epoch), th_est_idx] = ml_optimizer(th_est(p,k,epoch-1),optimizer,...
                    responses,trial_idx,p_correct,p_incorrect);
            else
                [th_est(p,k,epoch), th_est_idx] = ml_optimizer(th_est(p,k,epoch),optimizer,...
                    responses,trial_idx,p_correct,p_incorrect);
            end
            if isnan(th_est(p,k,epoch)) | isempty(th_est_idx)
                keyboard
            end
            optimizer_history(p,k,epoch) = optimizer;
            epoch = epoch+1;
        end
        response_accuracy(p,k) = mean(responses);
        emotion_history(p,k,:) = emo_strat;
        trial_idx_selected(p,k,:) = trial_idx;
        track110(p,k,:) = track_selected;
    end
end
%% Evaluate simulation results
for i = 1:test_length
    th_model = th_est(:,:,i);
    rho(i) = corr(th_model(:),th(:),'rows','pairwise');
end

figure
hold on
plot(rho,'LineWidth',5)
ylabel('Correlation Coefficient','FontSize',32);
xlabel('Test Length','FontSize',24);
set(gca,'FontSize',32,'LineWidth',2)
xlim([1 size(rho,2)])
xline(20,'LineWidth',5)
title('Correlation of true ability and test estimate')
box on
grid on
hold off

figure
hold on
scatter(th_est(1,:,20),th(1,:),100,'filled');
xlabel('Optimizer estimates','FontSize',24);
ylabel('True estimates','FontSize',24);
title({['Scatterplot of 1 permutation, TestLength=20'], ['r=' num2str(round(corr(th_est(1,:,20),...
    th(1,:)),2))]},'FontSize',28)
set(gca,'FontSize',32,'LineWidth',2)
h = lsline();
h.LineWidth = 5;h.Color = 'k';
xlim([-3 3]);
ylim([-3 3])
box on
grid on
hold off