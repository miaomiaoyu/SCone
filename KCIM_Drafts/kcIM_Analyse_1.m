% kcIM_Analyse2018

% Breaks up data into 1 second chunks and gets average power.
% Looks at whether the greatest power corresponds to
% You have to be in the Koniocellular folder to start.

% These are with the new triggers for 5, 12 and 16 Hz...
% think about these frequencies a bit more...
% 5, 12, 16?

% 30/10/17, M.Y.
% last edited 30/10/17, M.Y.

clear; close all;

currDir=mfilename('fullpath');
[p,n,e]=fileparts(currDir);
cd(p);

% ===== REMEMBER TO CHANGE THIS! =====
%cd('/Users/miaomiaoyu/GoogleDrive/Matlab_Toolboxes/Projects/Koniocellular');
cd(p);

%% FOR NOW:
cd('/Users/miaomiaoyu/GoogleDrive/Matlab_Toolboxes/Projects/Koniocellular');

%% THE SCRIPT
goodChans=1:65;
chansToAnalyse=31:33;
nGoodChans=length(goodChans);

% * note: these are triggers for 6, 12, 16 Hz.

% Triggers:
% Lum 5 -------- 24
% Lum 12 ------- 59
% Lum 16 ------- 79
% S Cone 5 ----- 34
% S Cone 12 ---- 83
% S Cone 16 ---- 111
% R-G 5 -------- 54
% R-G 12 ------- 131
% R-G 16 ------- 175
% isi ---------- 4
% sec ---------- 1
% pause -------- 22
% unpause ------ 29

analyseFolder='KCIM_EEG';
condTriggers=[24, 59, 79, 34, 83, 111, 54, 131, 175];

% EEG file directory
exptPath=[pwd, '/'];
eegDir=[exptPath, analyseFolder]; % Folder containing individual subj folders.
eegDataPath=dir(eegDir);
isDataFile=[eegDataPath.isdir];
subFolders=eegDataPath(isDataFile);

resultsDir=('./KCIM_AnalyseResults');

for thisFolderIndex=1:length(subFolders)
    dirPath{thisFolderIndex}=strcat(eegDir, '/', subFolders(thisFolderIndex).name);
    % dirPath gives you names of full folder directories: '.', '..', 'S1', 'S2', etc.
end

otherTriggers=[1, 4, 22, 29]; % Don't think I use this though.

for thisFolderIndex=3:length(dirPath) % This loops through each participant.
    
    fListDat=dir(fullfile(dirPath{thisFolderIndex},'*.cnt'));
    fListTrg=dir(fullfile(dirPath{thisFolderIndex},'*.evt'));
    
    % We loop over these files and load them in one by one.
    
    for thisSubjIndex=1:length(fListDat)
        fileName=fullfile(dirPath,fListDat(thisSubjIndex).name);
        EEG=mmy_extractEEGData(fileName{thisFolderIndex}, goodChans, 0);
        
    end
    
    if ~isempty(EEG)
        size(EEG);
        disp(EEG);
    else
        disp('Error extracting EEG data');
    end
%     
%     condTriggers=[24, 59, 79, 34, 83, 111, 54, 131, 175]; 
%     
%     if subFolders(thisFolderIndex).name == 'GBE'
%         condTriggers=[9, 59, 89, 13, 83, 125, 21, 131, 197];
%     end
    
    resultsfName = [resultsDir, '/', subFolders(thisFolderIndex).name];
    
    
    % =============================
    %% *** Tidying up the data ***
    
    % Noise reduction.
    
    % 1) I'll find the seconds where blink electrode shows great activity,
    % and chop those out.
    
    % 2) Take out the pause Index.
    
    % Here, I'm trying to cut off the 'excess' EEG recorded after the last ISI
    % trigger (indicating that the experiment has ended), and any chunks of
    % EEG data recorded during 'pause' periods. Although I'm not entirely
    % sure this is needed.
    
    pauseIndex = find(EEG.eventList == 22);
    unpauseIndex = find(EEG.eventList == 29);
    
    % *** insert pause code here ***
    
    exptEndTime = EEG.timeList(end) ; % The total time(s) from start to the last trigger.
    exptEndIndex = exptEndTime * EEG.rate;
    
    blinkPnts=EEG.data(end,:);
    
    EEG.data=EEG.data(chansToAnalyse,:);
    dataPnts = EEG.data; %(1:exptEndIndex);
    dataPnts=mean(dataPnts);
    
    % I'm not sure how much this helps, since we'll only be picking out the
    % 10-12 seconds coming after each condition code anyway.
    
    %% *** Main bit - picking out the data needed for FFT ***
    
    % This might not be the best way, but it's what I've thought to do:
    % First, I'll find the 15 indices (there were 15 reps in the expt) of each condition
    % trigger, then chop up and FFT the 10 seconds of data that comes after
    % the condition code.
    
    % I'll then average this across all three ppts, then try to present it
    % in a 3x3 figure.
    
    for triggerNo = 1:length(condTriggers) % Each one of these triggers is a differnt conditition e.g. 6Hz, Lum
        
        % First, find the index of these triggers (should be nReps long)
        triggerIndex = find(EEG.eventList == condTriggers(triggerNo));
        
        %disp(condTriggers(triggerNo));
        
        % Prep for if triggers found is less than rep no.
        if length(triggerIndex) < 15 % ZK only had 11 reps rather than 15.
            fprintf('Less than 15 triggers found for %g. ', condTriggers(triggerNo));   
        elseif length(triggerIndex) > 15 % cap it at 15
            triggerIndex = triggerIndex(1:15);  
        end
        
        nReps=length(triggerIndex); % just in case it's less than 15.
        
        nBins = 10; %7
        
        for repNo = 1:nReps
            
            for binNo = 1:nBins
                
                secondIndex(repNo, binNo)= triggerIndex(repNo) + binNo-1;
                
                % secondIndex is an array where the first element of every
                % row is the index of where the trigger appears (i), with
                % the subsequent 10 boxes being i + 1:10. I'm assuming that
                % the following 10 triggers following the cond code would
                % all be my '1' triggers: each denoting a second.
                % secondIndex
                
            end
        end
        
        dataPntIndex = (EEG.timeList(secondIndex) * 1/1000) * EEG.rate;
        % this is the corresponding time when that trigger was found.
        
        %dataPntIndex = floor((EEG.timeList(secondIndex)) * EEG.rate);
        % this is the corresponding time when that trigger was found.
        %dataPntsPerTrial=zeros(150,EEG.rate);
        
        for i = 1:numel(dataPntIndex)
            dataPntsPerTrial(i,:) = dataPnts(:,dataPntIndex(i):(dataPntIndex(i))+(EEG.rate-1));
            %dataPntsPerTrial(i,:) = dataPnts(dataPntIndex(i):(dataPntIndex(i)+(EEG.rate-1)));
        end % this is your reshapedWave essentially
        
        % Compute the raw variance of each bin to see if there are
        % very noisy bits...
        binPower=std(dataPntsPerTrial,[],2);
        %  figure(101);
        %  hist(binPower,20);
        badBins=find(zscore(binPower)>3); % Get rid of things more than x std from the mean
        dataPntsPerTrial(badBins,:)=0;
        
        
        range = 2:70; % we're not interested in much beyond that.
        
        ftDataPnts = fft(dataPntsPerTrial,[],2)/EEG.rate; % this is per condition. so each participant
        % should have 9 of these. Make sure we take the FT across dim 2
        % so subs are going down...
        
        %noiseRemoved = mmy_noiseRemoval(ftDataPnts, range); % REMOVES
        %1/F NOISE.
        
        meanPowerSpect = squeeze(mean(abs(ftDataPnts))); % ftDataPnts % incoh averaging?
        %meanPowerSpect = squeeze(abs(mean(ftDataPnts)));
        %meanPowerSpect = abs(squeeze(mean(ftDataPnts))); % COH
        
    end
    
    condPowerSpect(triggerNo,:) = meanPowerSpect;
    
end

plotTitle = {'Lum 6', 'Lum 12', 'Lum 16',...
    'S Cone 6', 'S Cone 12', 'S Cone 16',...
    'R-G 6', 'R-G 12', 'R-G 16'};

figure()

for i = 1:length(condTriggers)
    subplot(3,3,i)
    bar(condPowerSpect(i,range));
    h=title(plotTitle{i});
    set(h,'Visible','on');
    
end

save(resultsfName, 'condPowerSpect');


return

%% Last bit: averaging the data across all the participants

resultsDirPath = dir(resultsDir);
isResultsFile = [resultsDirPath.isdir];
resultsSubFolders = resultsDirPath(isResultsFile);

for thisFolderIndex = 3:length(resultsDirPath)
    load(resultsDirPath(thisFolderIndex).name);
    
    allResults(:,:,thisFolderIndex)=condPowerSpect;
    
end

avgResults=mean(allResults,3);

figure(100)

for i = 1:length(condTriggers)
    subplot(3,3,i)
    bar(avgResults(i,range));
    h=title(plotTitle{i});
    set(h,'Visible','on');
end



