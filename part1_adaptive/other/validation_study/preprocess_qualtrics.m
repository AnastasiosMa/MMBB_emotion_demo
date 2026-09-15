%%Correct the qualtrics missing ids issue, compile the final qualtrics file
pilot_participants = [23:31,33,34];
jspsych_dir = 'data/jatos_results/';
jatos_metadata = 'data/jspsych_metadata.json';
duplicate_participants = {'60fe9a0daa398bbc0d610dc0'};
pilot_ids = {'662fc8f6aa3d3b206672e28d','5c181b8d905b250001be42ca','66a27a561b8ed92d82f3b97a',...
    '5ea1a2e7939362055f0c325a','5dc9c2c0b9fad36ddca37632','663a3ee80ab62f61353a10f0','562a00dac8ffc20012513fbe',...
    '55d22025cc2b18000c0b9d9c','6715521230e88bbd98ef83ba','62d6d87d91cf24a3c1c0f10c','5d49d17b3dad1f0001e2aba1'};

qualtrics_data = readtable('qualtrics/EC_PRILITP_ERI_EN_Mavrolampados Anastasios_December 14, 2024_09.09_take1.csv');

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
qualtrics_data = qualtrics_data(:,{'UserID','score_facial','score_vocal'});
qualtrics_data(any(ismissing(qualtrics_data),2), :) = [];

qualtrics2 = readtable('qualtrics/ERI_part2_combined.csv');

microproms_data = readtable('data/proms_times.csv');
for i = 1:length(duplicate_participants)
    if any(strcmpi(microproms_data.Var3,duplicate_participants{i}))
        microproms_data.Var3{find(strcmpi(microproms_data.Var3,duplicate_participants{i}),1)}=[duplicate_participants{i},'_dble'];
    end
end
for i = 1:height(microproms_data)
    m_idx(i)=length(microproms_data{i,'Var3'}{1})>10;
end
microproms_data = microproms_data(m_idx,:);

prolific_demographics = readtable('data/prolific_export.csv');
for i = 1:length(duplicate_participants)
    if any(strcmpi(prolific_demographics.ParticipantId,duplicate_participants{i}))
        prolific_demographics.ParticipantId{find(strcmpi(prolific_demographics.ParticipantId,...
            duplicate_participants{i}),1)}=[duplicate_participants{i},'_dble'];
    end
end
%% jspsych
jatos_meta = fileread(jatos_metadata);
meta = jsondecode(jatos_meta);
for k = 1:length(meta.data.studyResults)
    participant_starttime(k,:) = [meta.data.studyResults{k}.id,...
        meta.data.studyResults{k}.startDate];
end
filenames = dir(jspsych_dir);
filenames = filenames(4:end); %HARDCODED FOR PILOT
jspsych_table = table();
for i = 1:length(filenames)
    filename_number = strsplit(filenames(i).name,'_');
    filename_number = filename_number{3};
    filename = [jspsych_dir, filenames(i).name, '/comp-result_',filename_number,'/data.txt'];
    r = fileread(filename);
    if any(pilot_participants==str2num(filename_number))
        r = r(1:length(r)/2); %HARDCODED FOR PILOT
    end
    jspsych_data = jsondecode(r);
    counter = 0;
    jspsych_table{i,'jatosID'} = str2num(filename_number);
    if any(pilot_participants==str2num(filename_number))
        jspsych_table{i,'UserID'} = {pilot_ids{find(pilot_participants==str2num(filename_number))}};
    else
        try
            jspsych_table{i,'UserID'} = {jspsych_data.trials{1}.userID};
        catch
            jspsych_table{i,'UserID'} = {'00'};
        end
    end
end
for i = 1:height(jspsych_table) %drop invalid ids
    valid_id(i)=length(jspsych_table.UserID{i})>10;
end
for i = 1:length(duplicate_participants)
if any(strcmpi(jspsych_table.UserID,duplicate_participants{i}))
   jspsych_table.UserID{find(strcmpi(jspsych_table.UserID,duplicate_participants{i}),1)}=[duplicate_participants{i},'_dble']; 
end
end
jatos_data = table();
for i = 1:length(meta.data.studyResults)
    jatos_data{i,'jatosid'} = meta.data.studyResults{i}.id;
    jatos_data{i,'StartDate'} = datetime(meta.data.studyResults{i}.startDate/1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    enddate = datetime(meta.data.studyResults{i}.endDate/1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    if ~isempty(enddate-jatos_data{i,'StartDate'})
        jatos_data{i,'EndDate'} = enddate;
        jatos_data{i,'jsDuration'} = seconds(enddate-jatos_data{i,'StartDate'});
    else
        jatos_data{i,'EndDate'} = NaT;
        jatos_data{i,'jsDuration'} = NaN;
    end
end
%concatenate jatos and jspsych
for i = 1:height(jspsych_table)
    if find(jspsych_table.jatosID(i)==jatos_data.jatosid)
       jatos_sorted(i,:) = jatos_data(find(jspsych_table.jatosID(i)==jatos_data.jatosid),:); 
    end
end
jspsych_data = [jspsych_table,jatos_sorted];
%% Keep only approved ids
%ids that have completed the study
approved_ids=prolific_demographics.ParticipantId(find(strcmpi(prolific_demographics.Status,'APPROVED')),:);
approved_ids = setdiff(approved_ids,qualtrics_data.UserID);
for i =1:length(approved_ids)
    approved_ids_prol(i,:) = prolific_demographics(find(strcmpi(prolific_demographics.ParticipantId,approved_ids{i})),:);
    try
        approved_ids_micro(i,:) = microproms_data(find(strcmpi(microproms_data.Var3,approved_ids{i})),:);
        approved_ids_js(i,:) = jspsych_data(find(strcmpi(jspsych_data.UserID,approved_ids{i})),:);
    end
end
approved_ids_micro{:,'promsDuration'} =  seconds(approved_ids_micro.DateLastAction - approved_ids_micro.DateStarted);
data = [approved_ids_js,approved_ids_micro,approved_ids_prol];
data.qualDuration = data.TimeTaken-(data.jsDuration+data.promsDuration);
data.DateStarted = data.DateStarted - seconds(3600);
data.DateLastAction = data.DateLastAction - seconds(3600);
writetable(data,'qualtrics/manual_check.csv')
%% Time difference matrix of microproms and qualtrics
for i=1:height(approved_ids_prol)
    st=approved_ids_prol.CompletedAt{i};
    time_difference = 7*3600;
    dt = datetime([st(1:10),' ',st(12:21)],'InputFormat','yyyy-MM-dd HH:mm:ss.S');
    time_d(:,i)= abs(seconds(dt-(qualtrics2.EndDate+seconds(time_difference))));
end
figure
imagesc(time_d)
colorbar