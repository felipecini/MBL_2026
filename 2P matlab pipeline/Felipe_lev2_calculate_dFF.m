function Felipe_lev2_calculate_dFF(dirsel,overwrite)
% KATIE_LEV2_CALCULATE_DFF finds dF/F.  Identifies artifact cells, and
% saves in datadFF struct as "goodCells"

% Katie Ferguson, Yale University, 2016

global info
global outputDirCardin
inputCa     = 'lev1_align_caframes';
analysis = 'lev2_calculate_dFF';
nDirs = length(info);

for iDir = dirsel
    exptag = info(iDir).dir;
    parts = split(exptag, '_');
    mouse = parts{1};
    outputDir = fullfile(outputDirCardin, analysis, mouse, exptag);
    mkdir(outputDir);
    
    dataDirCa       = fullfile(outputDirCardin, inputCa, mouse, exptag);
    dirInfoCa      = dir(dataDirCa);
    fileNamesCa    = {}; fileNamesState    = {}; fileNamesVis = {}; names = {};
    cnt = 0;
    nFiles     = length(dirInfoCa);
    
    
    for iFile = 1:nFiles
        if isempty(strfind(dirInfoCa(iFile).name, '.mat')), continue,end
        cnt = cnt + 1;
        names{cnt}     = dirInfoCa(iFile).name;
        fileNamesCa{cnt} = fullfile(dataDirCa, dirInfoCa(iFile).name);
    end
    
    for iFile = 1:length(fileNamesCa)
        
        outputFilename = fullfile(outputDir,names{iFile}(1:end-4));
        if exist([outputFilename '.mat'], 'file') && overwrite==0
            fprintf('skipping %s file %s \n', info(iDir).dir,outputFilename), continue,end
        
        fprintf('Processing %s\n', outputFilename);
        
        %load the Ca data
        load(fileNamesCa{iFile});
        
        %get rid of artifacts based on some rate of change threshold
        roc_artifact_thresh=300;
        rocCells=max(diff(datadFF.Ca,1,2),[],2)<roc_artifact_thresh;
        
        % define baseline as average of 10% min Ca points in a 10000 frame sliding window
        
        baselineProp=0.1; %percent of sorted signal to take for baseline
        F=datadFF.Ca;
        Fabsmin=min(min(F,[],1));
        if Fabsmin<0
            F=bsxfun(@plus,F,-1*Fabsmin);
        end
        
        windowsz=3000;    % window size setting, 100 second with framerate approximately 30Hz
        F0 = zeros(size(F));
        tic
        for idx=1:length(F)
            if idx <= windowsz/2
                block = F(:,1:idx+floor(windowsz/2));
            elseif idx > length(F)-floor(windowsz/2)
                block = F(:,idx-floor(windowsz/2):length(F));
            else
                block = F(:,idx-floor(windowsz/2):idx+floor(windowsz/2));
            end
            
            sortb = sort(block,2);
            F0(:,idx) = mean(sortb(:,1:floor(length(block)*baselineProp)),2);            % mean of lowest 10% of the 100 seconds sliding window
        toc
        end
        
        
        
        Fnum=bsxfun(@minus,F,F0);
        dFF=bsxfun(@rdivide,Fnum,F0);  %dFF=(F-F0)/F0
        
        datadFF.dFF=dFF;
        datadFF.goodCells=rocCells;
        
        save(outputFilename, 'datadFF');
        
    end
end








