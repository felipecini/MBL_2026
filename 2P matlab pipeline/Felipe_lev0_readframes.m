function Katie_lev0_readframes(dirsel, overwrite)
% KATIE_LEV0_READFRAMES reads in ca frames from spk2, and outputs struct called
% dataFrame with timeseries data, time of frame, number of Ca videos, 
% frame time of Ca video start/stop

% Katie Ferguson, Yale University, 2016

% 1/14/20 Calvin edited code to allow for smrx files to be read using the
% CED library


global outputDirCardin
global spike2Dir
global info

%load CED lib
cedpath='C:\Users\Peng\Desktop\meso_demo\MesoProcessing-master\process_spike2\CEDS64ML'; %location of cedpath
CEDS64LoadLib(cedpath);

% the analysis type, string to find data directories and name files
analysis  = 'lev0_readframes';

% loop over the various directories and read in the spike files
nDirs = length(info);
if nargin<1
    dirsel = 1:nDirs;
end
for iDir = dirsel
    exptag = info(iDir).dir;
    mouse = exptag(1:6);
    outputDir = fullfile(outputDirCardin, analysis, mouse, exptag);
    mkdir(outputDir);
    
    dataDir       = fullfile(spike2Dir, mouse, exptag);
    dirInfo      = dir(dataDir);
    fileNames    = {}; names = {};
    cnt = 0;
    nFiles     = length(dirInfo);
    for iFile = 1:nFiles
        if isempty(strfind(dirInfo(iFile).name, '.smr')), continue,end  %spk2 file
        cnt = cnt + 1;
        names{cnt}     = dirInfo(iFile).name;
        fileNames{cnt} = fullfile(dataDir, dirInfo(iFile).name);
    end
    
    dataAll = [];
    for iFile = 1:length(fileNames)
        fprintf('Processing %s\n', fileNames{iFile});
        outputFilename = fullfile(outputDir,names{iFile}(1:end-4));
        if exist([outputFilename 'mat'], 'file') && overwrite==0,
            fprintf('skipping %s file %s \n', info(iDir).dir,outputFilename), continue,end
        
        %begin Calvin's code using CED lib
        fhand = CEDS64Open(fileNames{iFile});
        %ichannum = min(CEDS64MaxChan(fhand),4);
        ichan = str2double(erase(info(iDir).frame,'chan')); %convert channel string from iDir into a number
        maxTimeTicks = CEDS64ChanMaxTime(fhand,1);
        dsec=CEDS64TimeBase(fhand); % sec per time tick, don't know if need
        [fRead,fVals,fTime] = CEDS64ReadWaveF(fhand,ichan,maxTimeTicks,0,maxTimeTicks); %don't know if need fRead and fTime yet
        data=fVals';
        time=linspace(0,length(data)/info(iDir).fsample,length(data)); %copied from Katie's code to calculate times
        %end Calvin's code using CED lib
        
        % find start and end indx and times
        %data_thresh=1000; 
        data_thresh = max(data)/2.25; %Really important, can result in +-2 frames different, need to pay attention 
        data_thresh1= max(data)/5;%added by Stefan, 08/25/2021, to avoid wrong number of frame, not sure if is suitable for all situation
        data_thresh2 = 1; %added by Stefan, 08/25/2021, to avoid wrong number of frame, not sure if is suitable for all situation
        gap_btw_movies=0.1 ; %min 100ms btw movies
        
       % frameInd=[intersect(find(data<data_thresh),find(diff(data)>data_thresh))]+1;
        frameInd=[intersect(find(data<data_thresh2),find(diff(data)>data_thresh1))]+1;%Edited by Stefan, 08/25/2021
        frameInd(find(diff(frameInd)==1)+1)=[]; %added by Stefan, 09/14/2021, to avoid consecutive frames.
        frameIndStart=[frameInd(1) frameInd(find(diff(frameInd)>(gap_btw_movies*info(iDir).fsample))+1)];
        frameIndEnd=[frameInd(find(diff(frameInd)>(gap_btw_movies*info(iDir).fsample))) frameInd(end)];
        
        
        
        %initialize dataFrame
        dataFrame=[];
        
        dataFrame.trial{1}=data;
        dataFrame.fsample=info(iDir).fsample;
        dataFrame.time{1}=time;
        dataFrame.label={names{iFile}(1:end-4)};
        dataFrame.cfg=[];
        dataFrame.time_timestamp{1}=dataFrame.time{1};
        dataFrame.frameInd=frameInd;
        dataFrame.imageInd(:,1)=frameIndStart;
        dataFrame.imageInd(:,2)=frameIndEnd;
                
        %find endpoint frames for videos in terms of low freq frame rate
        [~,frameEndpoints]=intersect(dataFrame.frameInd,dataFrame.imageInd);
        dataFrame.imageInd(:,3:4)=reshape(frameEndpoints,2,length(frameEndpoints)/2)';
        
        numMovies=length(frameIndStart);
        
        %checks if collect continuously.  Assumes either continuously or
        %discretely, not combination
        if numMovies==1 && max(dataFrame.imageInd(:,4)-dataFrame.imageInd(:,3))> info(iDir).framesPerMovie
            numMovies=ceil((dataFrame.imageInd(4)-dataFrame.imageInd(3))./info(iDir).framesPerMovie);
            dataFrame.imageInd=[]; %will re-assign
            movieInd=1:numMovies;
            dataFrame.imageInd(:,1)=frameInd(1+(movieInd-1).*info(iDir).framesPerMovie);
            movieIndEnd=movieInd.*info(iDir).framesPerMovie; 
            movieIndEnd(end)=size(frameInd,2); 
            dataFrame.imageInd(:,2)=frameInd(movieIndEnd);
            %find endpoint frames for videos in terms of low freq frame rate
            [~,frameEndpoints]=intersect(dataFrame.frameInd,dataFrame.imageInd);
            dataFrame.imageInd(:,3:4)=reshape(frameEndpoints,2,length(frameEndpoints)/2)';            
        end
        
        dataFrame.imageTime(:,1)=dataFrame.time{1}(frameIndStart);
        dataFrame.imageTime(:,2)=dataFrame.time{1}(frameIndEnd);
        
        dataFrame.imageIndLabel={'Start Index High Freq','End Index High Freq', 'Start Index Low Freq','End Index Low Freq'};
        dataFrame.imageTimeLabel={'Start Time','End Time'};
        dataFrame.numMovies=numMovies;
        % save the data to the disk
        save([outputFilename,'mat'], 'dataFrame','-v7.3');
    end
    
end
CEDS64CloseAll();
end




