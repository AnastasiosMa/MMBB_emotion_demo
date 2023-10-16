%% Export data matrices for CAT implementation
trial_info = readtable('../data/output/binary_responses/final_trial_info.csv');
rasch_mirt = readtable('../data/output/binary_responses/irt_models/final_rasch_mirt.csv');
mean_scores = readtable('../data/input/mean_ratings_set2.xls');
emo_idxs = [5,6,7,8,9];
emoNames = {'Angry','Fearful','Happy','Sad','Tender'};
item_difficulty = rasch_mirt{:,2};
item_emo = trial_info{:,2};
track = trial_info{:,end};
trialN = length(item_difficulty);
starting_item_difficulty = [1 -1];
%% export item table
example(1) = prctile(item_difficulty,25);
example(2) = prctile(item_difficulty,45);

idx_1 = find(item_emo==5);
idx_2 = find(item_emo==4);
%find example trials
[~,min_dist] = min(abs(example(1) - item_difficulty(idx_1)));
example_idx(1) = idx_1(min_dist);
[~,min_dist] = min(abs(example(2) - item_difficulty(idx_2)));
example_idx(2) = idx_2(min_dist);
example_data = trial_info(example_idx,:);

%remove example trials from data
trial_info(example_idx,:) = []; 
rasch_mirt(example_idx,:)  = [];
item_difficulty = rasch_mirt{:,2};
item_emo = trial_info{:,2};
track = trial_info{:,end};
trialN = length(item_difficulty);

%construct trial table
for i = 1:trialN
target_emo(i) = mean_scores{track(i),item_emo(i)};
soundtrack(i) = mean_scores{track(i),end-1};
end
trial_info.Properties.VariableNames = {'Label','Label1','Label2','Trials','Track110'};
trial_info.Difficulty = item_difficulty;
trial_info.TargetEmo= target_emo';
trial_info.Soundtrack= soundtrack';

%% export data matrices
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
%% export all
json_p_star = jsonencode(p_star');
json_information = jsonencode(information_test');

%save pstar
fid = fopen('../data/output/binary_responses/cat_data/pstar.json','w');
fprintf(fid,'%s',json_p_star);
fclose(fid);

%save information
fid = fopen('../data/output/binary_responses/cat_data/information.json','w');
fprintf(fid,'%s',json_information);
fclose(fid);

Xc = rows2vars(trial_info);
Xc = Xc(:,2:end);
json_trial_info = jsonencode(trial_info);
for i = 1:size(trial_info,2)
    t{i} = trial_info{:,i};
end
json_t = jsonencode(t);
%save trial info
fid = fopen('../data/output/binary_responses/cat_data/trial_info.json','w');
fprintf(fid,'%s',json_t);
fclose(fid);
%save labels
labels = trial_info.Properties.VariableNames;
writetable(array2table(labels),'../data/output/binary_responses/cat_data/trial_info_labels.txt')