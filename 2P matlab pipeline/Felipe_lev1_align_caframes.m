function Katie_lev1_align_caframes(dirsel, overwrite)
% KATIE_LEV1_ALIGN_CAFRAMES reads in Ca and frame data, and aligns appropriately.  

% Katie Ferguson, Yale University, 2016
% Modified by Felipe Cini, Yale University, 2026

global outputDirCardin
global info

% the analysis type, string to find data directories and name files
inputCa  = 'lev0_readca';
inputFrames  = 'lev0_readframes';
analysis = 'lev1_align_caframes';

% loop over the various directories and read in the ca files
nDirs = length(info);
if nargin<1
  dirsel = 1:nDirs;
end

for iDir = dirsel
    exptag = info(iDir).dir;
    parts = split(exptag, '_');
    mouse = parts{1};
    outputDir = fullfile(outputDirCardin, analysis, mouse, exptag);
    mkdir(outputDir);
    
    dataDirFrames       = fullfile(outputDirCardin, inputFrames, mouse, exptag);
    dirInfoFrames      = dir(dataDirFrames);
    fileNamesFrames    = {}; namesFrames = {};
    nFilesFrames     = length(dirInfoFrames);
    
    dataDirCa       = fullfile(outputDirCardin, inputCa, mouse, exptag);
    dirInfoCa      = dir(dataDirCa);
    nFilesCa     = length(dirInfoCa);
    CafileName  = {dirInfoCa(:).name};
    
    for iFileFrames= 1:nFilesFrames
        if isempty(strfind(dirInfoFrames(iFileFrames).name, '.mat')), continue,end  
        namesFrames    = dirInfoFrames(iFileFrames).name;
        fileNamesFrames= fullfile(dataDirFrames, dirInfoFrames(iFileFrames).name);
    end
    
    datadFF = [];
    
    fprintf('Processing %s\n', fileNamesFrames);
    outputFilename = fullfile(outputDir,namesFrames(1:end-4));
    if exist([outputFilename '.mat'], 'file') && overwrite==0
        fprintf('skipping %s file %s \n', info(iDir).dir,outputFilename), continue,end
    
    %load the frame data
    load(fileNamesFrames,'dataFrame');
    indx   = strfind(namesFrames,'_');
    sessionNum = str2num(namesFrames(indx(end)+1:indx(end)+5));  %%will need to change for Higley lab
    
    %load Ca data
    fileNamesCa    = {}; namesCa = {};
    dFF_Ind=[]; dFF_Time=[]; dFF_data=[];
    dFF_StartInd_hsamp=[]; dFF_EndInd_hsamp=[]; dFF_StartTime_hsamp=[]; dFF_EndTime_hsamp=[];
    dFF_Ind=[]; dFF_StartInd=[]; dFF_EndInd = [];
    
    iFileCa=1;
    namesCaCheck=namesFrames(1:indx(end));
    for iallCa=1:nFilesCa
        if isempty(strfind(CafileName{iallCa}, namesCaCheck)); continue,end  %spk2 file
        namesCa    = CafileName{iallCa};
        fileNamesCa= fullfile(dataDirCa, CafileName{iallCa});
    end
    
    if isempty(fileNamesCa)
        fprintf('No Calcium file for movie %d \n'); continue,
    end
    
    load(fileNamesCa,'dataCa');
    
    CaStartInd = dataFrame.imageInd(iFileCa,3);
    numFrames=dataFrame.imageInd(end,4);
    if numFrames<length(dataCa.trial{1})
        fprintf('Warning: %d fewer frames than Calcium points \n',length(dataCa.trial{1})-numFrames)
    elseif numFrames>length(dataCa.trial{1})
        fprintf('Warning: %d more frames than Calcium points \n',numFrames-length(dataCa.trial{1}))
    end
    
    dFF_StartInd(iFileCa) =CaStartInd + dataCa.time_start_ind{1} -1;
    dFF_EndInd(iFileCa) = dFF_StartInd(iFileCa) + numFrames -1 ;
    
    dFF_Ind=[dFF_Ind dataFrame.frameInd(dFF_StartInd(iFileCa):dFF_EndInd(iFileCa))]; %ind of every dFF frame
    
    dFF_data=[dFF_data dataCa.trial{1}];
    
    
    if size(dFF_data,1)>0  %if not all cells had artifact..
        dFF_Time=dataFrame.time{1}(dFF_Ind);  %time of every dFF frame
        
        datadFF.frameTime = dataFrame.time{1};
        datadFF.frameInd = dataFrame.frameInd;
        datadFF.caTime = dataCa.time{1};
        datadFF.time = dFF_Time;
        datadFF.ind = dFF_Ind;
        datadFF.Ca = dFF_data;
        datadFF.spks = dataCa.spks;
        datadFF.iscell = dataCa.iscell;
        datadFF.frameIntervalInd = dataFrame.imageInd;
        datadFF.frameIntervalTime = dataFrame.imageTime;
        datadFF.frameIntervalInd_label = {'Start Ind High Samp','End Ind High Samp','Start Frame Ind','End Frame Ind'};
        datadFF.frameIntervalTime_label = {'Start Frame Time','End Frame Time'};
        
        imageLength=datadFF.frameIntervalInd(:,4)-datadFF.frameIntervalInd(:,3);
        datadFF.imageInd=[[1;cumsum(imageLength(1:end-1)+1)], cumsum(imageLength+1)+1];
        datadFF.imageInd_label={'Start Image Ind','End Image Ind'};  %wrt Ca time only
        
        % save the data to the disk
        save(outputFilename, 'datadFF','-v7.3');
    end
end
    

  
  
  
  


