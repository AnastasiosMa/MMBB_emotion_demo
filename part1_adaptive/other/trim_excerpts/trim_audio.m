%trim audio files to the last 4 seconds and add 1 sec fade in
addpath(genpath('~/Documents/extra_Matlab_functions'))
cd mp3
dirNames = dir;
dirNames = dirNames(4:end);
sr = 44100;
stimuliThres = 10;
for i = 1:length(dirNames)
    a = mirgetdata(miraudio(dirNames(i).name,'Mono',0));
    a = squeeze(a);
    idx = 1:sr*stimuliThres;
    idx = int64(idx);
    trimmedAudioData{i} = a(idx,:);
    filter = linspace(1,0,sr*1);
    filter = [ones(1,sr*[stimuliThres-1]),filter];
    trimmedAudioData{i} = trimmedAudioData{i}.*[filter',filter'];
    dirNames(i).name = [dirNames(i).name(1:3) '.wav'];
end
cd ..
cd trimmed
for i = 1:length(dirNames)
    audiowrite(dirNames(i).name,trimmedAudioData{i},sr)
    %mp3write(trimmedAudioData{i},sr,16,dirNames(i).name,2)
end