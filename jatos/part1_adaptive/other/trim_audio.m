%trim audio files to the last 4 seconds and add 1 sec fade in
addpath(genpath('~/Documents/extra_Matlab_functions'))
cd part2audio/
dirNames = dir;
dirNames = dirNames(4:end);
sr = 44100;
stimuliThres = 8;
stimuliCut = 'start';
for i = 1:length(dirNames)
    a = mirgetdata(miraudio(dirNames(i).name,'Mono',0));
    a = squeeze(a);
    if strcmpi(stimuliCut,'start')
        idx = 1:sr*stimuliThres;
    else
        idx = length(a)-sr*stimuliThres:length(a);
    end
    idx = int64(idx);
    trimmedAudioData{i} = a(idx,:);
    if strcmpi(stimuliCut,'start')
        filter = linspace(1,0,sr*1);
        filter = [ones(1,sr*[stimuliThres-1]),filter];
    else
        filter = linspace(0,1,sr*1);
        filter = [filter,ones(1,sr*[stimuliThres-1])];
    end
    trimmedAudioData{i} = trimmedAudioData{i}.*[filter',filter'];
    dirNames(i).name = [dirNames(i).name(1:3) '.wav'];
end
cd ..
cd trimmed_part2_wav
for i = 1:length(dirNames)
    audiowrite(dirNames(i).name,trimmedAudioData{i},sr)
    %mp3write(trimmedAudioData{i},sr,16,dirNames(i).name,2)
end