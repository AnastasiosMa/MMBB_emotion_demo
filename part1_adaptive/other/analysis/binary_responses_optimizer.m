data = readtable('../data/output/binary_responses/fear_binary_responses_only.csv');
trial_info = readtable('../data/output/binary_responses/fear_trial_info.csv');
trialN = size(data,2);
emoNames = {'Angry','Fearful','Happy','Sad','Tender'};
difficulty_item = mean(data{:,1:trialN},'omitnan');
rasch_mirt = readtable('../data/output/binary_responses/irt_models/rasch_mirt.csv');
%% Optimizer denominator
%create probability of correct and wrong sample answers
theta_step = 0.05;
theta_low = -8;
theta_high = 8;
thetaRange = theta_low:theta_step:theta_high;
for j = 1:trialN %trial
    k=1;
    for th = thetaRange
        p_correct(j,k) = 1/(1+exp(-(th-rasch_mirt{j,2})));
        p_incorrect(j,k) = 1-(1/(1+exp(-(th-rasch_mirt{j,2}))));
        k = k+1;
    end
end
%i = th_idx, j = item order, u = participant responses th = theta value
%index optimizer formula
optimizer = @(i,u,j) sum(u-p_correct(j,i)')/(-sum(p_correct(j,i)'.*p_incorrect(j,i)'));
%% Test optimizer
init_th = 0; %initial starting point of theta
learning_rate = 0.5;
min_criterion = theta_step; %optimizer stops if Δθ < criterion
n_iter = 500; %number of iterations before manual stop
%timings = nan(size(data,1),size(data,2));
iterations = nan(size(data,2),n_iter);
optimizer_training_history = cell(1,size(data,1));
tic
for participant = 1:size(data,1) %per participant
    optimizer_training_history{participant} = iterations;
    participant_data = data{participant,:};
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
            elseif th>theta_high
                th_idx = length(thetaRange);
            end
            %find discrete index for theta within range
            th_idx = find(th<=thetaRange+theta_step & th>thetaRange);
            if isempty(th_idx)
                %keyboard
            end
            delta_th = optimizer(th_idx,participant_data(p_resp(1:trial)),p_resp(1:trial));
            th = th-delta_th;
            optimizer_training_history{participant}(trial-1,iter) = th;
            iter = iter+1;
        end
        %timings(participant,trial-1) = toc(tStart);
    end
end
toc