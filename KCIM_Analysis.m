% kcIM_Analysis

% Breaks up data into 1 second chunks and gets average power.
% Looks at whether the greatest power corresponds to
% You have to be in the Koniocellular folder to start.

% These are with the new triggers for 5, 12 and 16 Hz...
% think about these frequencies a bit more...
% 5, 12, 16?

% 30/10/17, M.Y.
% last edited: 25/01/17, M.Y.

%% General Directories (This bit never changes):

clear; close all;

thisComputer = computer;

if strcmp(thisComputer, 'MACI64')
    curDir = ('/Users/miaomiaoyu/Documents/GitHub/NeuralOscillations');
else
    curDir = ('/wadelab_shared/Projects/NeuralOscillations');
end

%% Additional Directories (..:: This is the bit you change! ::..)

projectFolder = '/Koniocellular'; % KCIM, Rodent, or Psychophysics
dataFolder = '/KCIM_EEG'; % Folder that has S1, S2, S3 subfolders.
resultsFolder = '/KCIM_Analysis_Results';

%% General Directories

curDir = [curDir, projectFolder];

dataDir =  [curDir, dataFolder];
resultsDir = [curDir, resultsFolder];
dataPath = dir(dataDir);
isDataFile = [dataPath.isdir];
subFolders = dataPath(isDataFile);

% dirPath gives you full folder directories ('.', '..', 'S1', 'S2' etc.)
for thisFolderIndex = 1:length(subFolders)
    dirPath{thisFolderIndex} = strcat(dataDir, '/', subFolders(thisFolderIndex).name);
end

goodChans=1:65;
chansToAnalyse=31:33;
nGoodChans=length(goodChans);

for thisFolderIndex=3:length(dirPath) % This loops through each participant.
    
    fListDat=dir(fullfile(dirPath{thisFolderIndex},'*.cnt'));
    fListTrg=dir(fullfile(dirPath{thisFolderIndex},'*.evt'));
    
    % We loop over these files and load them in one by one.
    
    for thisSubjIndex=1:length(fListDat)
        fileName=fullfile(dirPath,fListDat(thisSubjIndex).name);
        EEG=mmy_Extract_Subject_Data_EEG_ANT(fileName{thisFolderIndex}, goodChans, 0);
        
    end
    
    if ~isempty(EEG)
        size(EEG);
        disp(EEG);
    else
        disp('Error extracting EEG data');
    end
    
    
    %% THE SCRIPT:
    % * note: these are triggers for 5, 12, 16 Hz.
    
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
    
    condTriggers=[24, 59, 79, 34, 83, 111, 54, 131, 175];
    
    if strcmpi(subFolders(thisFolderIndex).name, 'GBE')
        condTriggers=[9, 59, 89, 13, 83, 125, 21, 131, 197];
    end    % GBE has a different set of triggers.
    
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
    % You should pick out the bits where the pause code was used and
    % discard the data between... find pause code, and the next unpause
    % code, discard the data, or just discard the entire rep it comes with.
    % find big blinkPnts and discard those as well.
    
    nWholeSeconds = floor(EEG.timeList(end) * 1/1000); % The total time (s) from start to the last trigger.
    
    nDataPoints = (nWholeSeconds+1) * EEG.rate; % add 1 second, just in case sometimes the last second run a little further.
    
    %% Take out the blink noises!
    
    zsLimit = 3.8; % Looks most suitable according to figures.
    displayFig = 0; % Makes mmy_Noise_Extraction_Zscore enormously slow.
    binSize = 400; % The average blink is about 400ms...
    
    nRemData = rem(nDataPoints, binSize);
    nDataPoints = nDataPoints - nRemData;
    
    dataPoints=mean(EEG.data(chansToAnalyse, :));
    floorDataPoints=mean(EEG.data(chansToAnalyse, 1:nDataPoints));
    
    blinkPoints=EEG.data(end, 1:nDataPoints);
    
    %zsIndices = mmy_Noise_Extraction_Zscore(blinkPoints, zsLimit, displayFig);
    %sumBps = mmy_Noise_Per_Bin(zsIndices, zsLimit, nDataPoints, EEG.rate/2, exptEndTime*2);
    
    [trashBits, zsIndices, sumBPS] = mmy_Noise_Extraction_Zscore(blinkPoints, zsLimit, nDataPoints, ...
    binSize, nDataPoints/binSize, displayFig);
    
    floorDataPoints(trashBits) = nan;
    
    % If the error's regarding reshape elements must not change - check
    % that you've got the right binSize (y/2) and exptEndTime (z*2).
    
    
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
        
        % Prep for if triggers found is less than rep no.
        
        if length(triggerIndex) < 15 % ZK only had 11 reps rather than 15.
            fprintf('Less than 15 triggers found for %g. ', condTriggers(triggerNo));
        elseif length(triggerIndex) > 15 % cap it at 15
            triggerIndex = triggerIndex(1:15);
        end
        
        nReps=length(triggerIndex); % just in case it's less than 15.
        
        nBins = 8; % It's only really 8 bins we have for some reason - so stick with this.
        
        for repNo = 1:nReps
            
            for binNo = 1:nBins
                
                secondIndex(repNo, binNo)= (triggerIndex(repNo)) + binNo-1;
                
                % secondIndex is an array where the first element of every
                % row is the index of where the trigger appears (i), with
                % the subsequent 10 boxes being i + 1:10. I'm assuming that
                % the following 10 triggers following the cond code would
                % all be my '1' triggers: each denoting a second.
                % secondIndex
                
                % in reality, there are only eight '1's that follow my main
                % condition trigger... also the last 1 coincides with the
                % isi 4 in terms of timing - so perhaps we only have seven.
                
                
            end
            
        end
        
        for repNo = 1:nReps
            
            dataPntIndex = round(EEG.timeList(secondIndex) * EEG.rate/1000);
            % this is the corresponding time when that trigger was found.
            
            % dataPntIndex = floor((EEG.timeList(secondIndex)) * EEG.rate);
            
            for i = 1:numel(dataPntIndex)
                dataPntsPerTrial(i,:) = floorDataPoints(dataPntIndex(i):(dataPntIndex(i)+(EEG.rate-1)));
            end % this is your reshapedWave essentially
            
            % Compute the raw variance of each bin to see if there are
            % very noisy bits...
            
            binPower=std(dataPntsPerTrial,[],2);
            %  figure(101);
            %  hist(binPower,20);
            %  badBins=find(zscore(binPower)>3); % Get rid of things more than x std from the mean
            %  dataPntsPerTrial(badBins,:)=0;
            
            myrange = 2:70; % we're not interested in much beyond that.
            
            ftDataPnts = fft(dataPntsPerTrial,[],2)/EEG.rate;
            % this is per condition. so each participant
            % should have 9 of these. Make sure we take the FT across dim 2
            % so subs are going down...
            
            incohPowerSpect = squeeze(nanmean(abs(ftDataPnts))); % incoh averaging
            cohPowerSpect = abs(squeeze(nanmean(ftDataPnts))); % coh averaging
            
        end
        
        allIncohPowerSpect(triggerNo,:) = incohPowerSpect;
        allCohPowerSpect(triggerNo,:) = cohPowerSpect;
        
        % disp(condTriggers(triggerNo));
        
    end
    
    maxYLim = max(allCohPowerSpect(myrange), [], 1);
    minYLim = min(allCohPowerSpect(myrange), [], 1);
    
    plotTitle = {'Lum 5', 'Lum 12', 'Lum 16',...
        'S Cone 5', 'S Cone 12', 'S Cone 16',...
        'R-G 5', 'R-G 12', 'R-G 16'};
    
    displayFigure = 1;
    
    if displayFigure
        
        figure()
        
        for i = 1:length(condTriggers)
            subplot(3,3,i)
            bar(allIncohPowerSpect(i,myrange));
            h=title(plotTitle{i});
            set(h,'Visible','on');
            
        end
        
        figure()
        
        for i = 1:length(condTriggers)
            subplot(3,3,i)
            bar(allCohPowerSpect(i,myrange));
            h=title(plotTitle{i});
            set(h,'Visible','on');
            %ylim([minYLim, maxYLim ]);
            
        end
        
    end
    
    save([resultsfName, '_InCoh'], 'allIncohPowerSpect');
    save([resultsfName, '_Coh'], 'allCohPowerSpect');
    
end

% return
% - %

%% Last bit: averaging the data across all the participants

resultsDirPath = dir(resultsDir);
isResultsFile = ~[resultsDirPath.isdir];
resultsSubFolders = resultsDirPath(isResultsFile);

%clear('allCohResults'); clear('allIncohResults'); % when coding

cohFileIndex=1;
incohFileIndex=1;

for thisFolderIndex = 1:length(resultsSubFolders)
    
    load([resultsDir, '/', resultsSubFolders(thisFolderIndex).name]); % load in the data files 1 by 1
    
    [rPath, rName, rExt]=fileparts(resultsSubFolders(thisFolderIndex).name);
    
    if contains(rName, '_Coh')
        allCohResults(:,:,cohFileIndex) = allCohPowerSpect;
        %disp('coh')
        cohFileIndex = cohFileIndex+1;
        
    elseif contains(rName, '_InCoh')
        allIncohResults(:,:,incohFileIndex)=allIncohPowerSpect;
        %disp('incoh')
        incohFileIndex=incohFileIndex+1;
        
    end
end

avgCohResults=mean(allCohResults, 3);
avgIncohResults=mean(allIncohResults, 3);

figure(100)
for i = 1:length(condTriggers)
    subplot(3,3,i)
    bar(avgCohResults(i, myrange));
    h=title(plotTitle{i});
    set(h, 'Visible', 'on');
end

figure(101)
for i = 1:length(condTriggers)
    subplot(3,3,i)
    bar(avgIncohResults(i, myrange));
    h=title(plotTitle{i});
    set(h, 'Visible', 'on');
end




