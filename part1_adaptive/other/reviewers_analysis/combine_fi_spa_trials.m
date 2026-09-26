%% Concatenate the common finnish and spanish trials
binary_responses_fi = readtable('data/binary_responses_fi_all.csv');
trial_info_fi = readtable('data/trial_info_fi_all.csv');

binary_responses_spa = readtable('data/binary_responses_spa_all.csv');
trial_info_spa = readtable('data/trial_info_spa_all.csv');

trial_info_fi.Index = strcat(string(trial_info_fi.TargetEmo), string(trial_info_fi.ComparisonEmo),string(trial_info_fi.Track110));
trial_info_spa.Index = strcat(string(trial_info_spa.TargetEmo), string(trial_info_spa.ComparisonEmo),string(trial_info_spa.Track110));

commonIDs = intersect(trial_info_fi.Index, trial_info_spa.Index);

fi_keep = ismember(trial_info_fi.Index, commonIDs);
spa_keep = ismember(trial_info_spa.Index, commonIDs);

trial_info_fi = trial_info_fi(fi_keep, :);
trial_info_spa = trial_info_spa(spa_keep, :);

binary_responses_spa = binary_responses_spa(:,spa_keep);
binary_responses_fi = binary_responses_fi(:,fi_keep);

[trial_info_fi,idx] = sortrows(trial_info_fi,'Index'); 
binary_responses_fi = binary_responses_fi(:,idx);

[trial_info_spa,idx] = sortrows(trial_info_spa,'Index'); 
binary_responses_spa = binary_responses_spa(:,idx);

writetable(binary_responses_fi,'data/binary_responses_fi.csv');
writetable(trial_info_fi,'data/trial_info_fi.csv');

writetable(binary_responses_spa,'data/binary_responses_spa.csv');
writetable(trial_info_spa,'data/trial_info_spa.csv');