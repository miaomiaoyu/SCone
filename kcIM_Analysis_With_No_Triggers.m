% kcIM_Analysis_With_T_Tests

% Breaks up data into 1 second chunks and gets average power.
% Looks at whether the greatest power corresponds to
% You have to be in the Koniocellular folder to start.

% These are with the new triggers for 5, 12 and 16 Hz...
% think about these frequencies a bit more...
% 5, 12, 16?

% 30/10/17, M.Y.
% last edited: 25/01/17, M.Y.


% improvements to be made: use ICA as blink electrode noise removal
% technique

% run the data through the classification function too

%

%% General Directories (This bit never changes):

clear; close all;gramm

 computeName=char(java.net.InetAddress.getLocalHost.getHostName);
if strcmp(computeName,'d2') % Are we on D2?

    EEGpath = '/wadelab_shared/Projects/NeuralOscillations//';

elseif strcmp(computeName, 'nas10-240-125-18.york.ac.uk')  % Are we on Miaomiao's mac?
    EEGpath = ('/Users/miaomiaoyu/Documents/GitHub/NeuralOscillations');
% Assume we are at YNiC
else
    EEGpath = '/groups/labs/wadelab/data/Miaomiao/NeuralOscillations//';

end


% thisComputer = computer;
% 
% if strcmp(thisComputer, 'MACI64')
%     curDir = ('/Users/miaomiaoyu/Documents/GitHub/NeuralOscillations');
% else
%     curDir = ('/wadelab_shared/Projects/NeuralOscillations');
% end

curDir=EEGpath;

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
    %
    %     zsIndices = mmy_Noise_Extraction_Zscore(blinkPoints, zsLimit, displayFig);
    %     sumBps = mmy_Noise_Per_Bin(zsIndices, zsLimit, nDataPoints, EEG.rate/2, exptEndTime*2);
    %
    %     %[trashBits, zsIndices, sumBPS]=mmy_Noise_Extraction_Zscore(blinkPoints, zsLimit, nDataPoints, ...
    %     %binSize, nDataPoints/binSize, displayFig);
    %
    %      floorDataPoints(trashBits) = nan;
    
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
        
        nBins = 12; % It's only really 8 bins we have for some reason - so stick with this.
        
        %for repNo = 1:nReps
        
        % for binNo=1:nBins
        
        %secondIndex(repNo)= (triggerIndex(repNo));
        
        % secondIndex is an array where the first element of every
        % row is the index of where the trigger appears (i), with
        % the subsequent 10 boxes being i + 1:10. I'm assuming that
        % the following 10 triggers following the cond code would
        % all be my '1' triggers: each denoting a second.
        % secondIndex
        
        % in reality, there are only eight '1's that follow my main
        % condition trigger... also the last 1 coincides with the
        % isi 4 in terms of timing - so perhaps we only have seven.
        
        %  end
        
        
        %end
        
        for repNo = 1:nReps
            
            dataPntIndex = round(EEG.timeList(triggerIndex) * EEG.rate/1000);
            % these are the actual data point indices of the first twelve,
            % so essentially if dataPntIndex = 234, you want the next few
            % to be 1234, 2234, 3234... etc. it's dataPntIndex + 1000*i
            % (with i = [1:12]-1).
            
            % this is the corresponding time when that trigger was found.
            
            % dataPntIndex = floor((EEG.timeList(secondIndex)) * EEG.rate);
            
            
            %             for i = 1:numel(dataPntIndex)
            %                 dataPntsPerTrial(i,:) = floorDataPoints(dataPntIndex(i):(dataPntIndex(i)+(EEG.rate-1)));
            %             end % this is your reshapedWave essentially
            %
            
            for i = 1:length(dataPntIndex)
                for binNo=1:11
                    dataPntsIndexAll(i,binNo)=(dataPntIndex(i)+(EEG.rate*binNo)-EEG.rate);
                    
                    % this is 15 repetitions of 12 seconds long each
                    % time... (so all of 1 condition they would see for the
                    % whole experiment)...
                end
            end
            
            for i = 1:numel(dataPntsIndexAll)
                dataPntsPerTrial(i,:) = floorDataPoints(dataPntsIndexAll(i):(dataPntsIndexAll(i)+(EEG.rate-1)));
            end
            
            % Compute the raw variance of each bin to see if there are
            % very noisy bits...
            
            binPower=std(dataPntsPerTrial,[],2);
            %  figure(101);
            %  hist(binPower,20);
            %  badBins=find(zscore(binPower)>3); % Get rid of things more than x std from the mean
            %  dataPntsPerTrial(badBins,:)=0;
            
            myrange = 2:100; % we're not interested in much beyond that.
            
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
    
    
    plotTitle = {'Lum 5', 'S Cone 5', 'RG 5',...
        'Lum 12', 'S Cone 12', 'RG 12',...
        'Lum 16', 'S cone 16', 'RG 16'};
    
    
    displayFigure = 0;
    
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
    
    if strfind(rName, '_Coh')
        allCohResults(:,:,cohFileIndex) = allCohPowerSpect;
        cohFileIndex = cohFileIndex+1;
        
    elseif strfind(rName, '_InCoh')
        allIncohResults(:,:,incohFileIndex)=allIncohPowerSpect;
        incohFileIndex=incohFileIndex+1;
        
    end
end

% for subjNo = 1:(16-1)
%     diffCoh(:,:,subjNo)=allCohResults(:,:,subjNo)-allCohResults(:,:,subjNo+1);
% end

avgCohResults=mean(allCohResults, 3);
avgIncohResults=mean(allIncohResults, 3);

%%             Figures for Publishing:
%
%     |   LUM 5   |   SCONE 5   |   RG 5   |
%     |   LUM 12  |   SCONE 12  |   RG 12  |
%     |   LUM 16  |   SCONE 16  |   RG 16  |

plotTitle = {'Lum 5', 'S Cone 5', 'RG 5',...
    'Lum 12', 'S Cone 12', 'RG 12',...
    'Lum 16', 'S cone 16', 'RG 16'};

% for thisPT=1:9
%     newPT{thisPT}=plotTitle{orderForFig(thisPT)};
% end

orderForFig=[1 4 7 2 5 8 3 6 9];

avgCohResults=avgCohResults(orderForFig,:);
avgIncohResults=avgIncohResults(orderForFig,:);

figure(100)

for i = 1:length(condTriggers)
    s=subplot(3,3,i);
    bar(avgCohResults(i, myrange));
    h=title(plotTitle{i});
    set(h, 'Visible', 'on');
    
    linkaxes([s]);
    ylim([0 0.6])
    
end

%saveas(gcf, 'Fig_100_Coh_100Hz.pdf');

figure(101)

for i = 1:length(condTriggers)
    s=subplot(3,3,i);
    bar(avgIncohResults(i, myrange));
    h=title(plotTitle{i});
    set(h, 'Visible', 'on');
    linkaxes([s]);
    ylim([0 1])
end

%saveas(gcf, 'Fig_101_Incoh_100Hz.pdf');

%% THIS IS THE BIT YOU CHANGE!

figure(102);

outputRange=2:100; % used to be 2:150

rsCoh=reshape(avgCohResults, [3,3,1000]);
rsIncoh=reshape(avgIncohResults, [3,3,1000]);

% avgCohResults               %  rsCoh
% ------------------------------------------------
% lum 5                 lum 5  | scone 5  |  RG 5
% lum 12                lum 12 | scone 12 |  RG 12
% lum 16                lum 16 | scone 16 |  RG 16
% scone 5
% scone 12
% scone 16
% RG 5
% RG 12
% RG 16

compA=[3 2 3];   % i.e. 3-1 => RG - lum
compB=[1 1 2];   %      2-1 => S - lum
                 %      3-2 => RG - S

for freqNo=1:3 % we'll only deal with one frequency at a time
    for compNo=1:3  % subtracting powers from diff colours from each other
        diffPowerCoh(freqNo,compNo,:)=rsCoh(freqNo,compA(compNo),outputRange)-rsCoh(freqNo,compB(compNo),outputRange);
        diffPowerIncoh(freqNo,compNo,:)=rsIncoh(freqNo,compA(compNo),outputRange)-rsIncoh(freqNo,compB(compNo),outputRange);
    end
end
        
freqTitle={'6 Hz', '12 Hz', '18 Hz'};
diffTitle={'RG - lum', 'S - lum', 'RG - S'};

for freqNo=1:3
    figure()     % 1-3 Coh: produces 2 figures showing differences..
    h=title(freqTitle{freqNo});
    set(h, 'Visible', 'on');
    for compNo=1:3
        subplot(3,1,compNo)
        bar(squeeze(diffPowerCoh(freqNo,compNo,:)));
        h=title(diffTitle{compNo});
        set(h, 'Visible', 'on');
    end
end

for freqNo=1:3
    figure()     % 4-6 Incoh
    h=title(freqTitle{freqNo});
    set(h, 'Visible', 'on');
    for compNo=1:3
        subplot(3,1,compNo)
        bar(squeeze(diffPowerIncoh(freqNo,compNo,:)));
        h=title(diffTitle{compNo});
        set(h, 'Visible', 'on');
    end
end

sumPowerCoh=squeeze(sum(avgCohResults(:,outputRange),2));
sumPowerIncoh=squeeze(sum(avgIncohResults(:,outputRange),2));

figure(105);

subplot(2,1,1)
im=imagesc(reshape(sumPowerCoh,[3,3]));
colormap(flipud(gray(256)));
colorbar;

%j= insertText(im, [100 315], 'pepe');
%imshow(j);

subplot(2,1,2)
im=imagesc(reshape(sumPowerIncoh,[3,3]));

colormap(flipud(gray(256)));
colorbar;

%% Grouping powers by endogenous frequency range

inputFreq=[5, 12, 16];
maxFreq = 150; % the maximum freq we're working with here.

reshapeCoh=reshape(allCohResults, [3,3,1000,16]);
reshapeIncoh=reshape(allIncohResults, [3,3,1000,16]);

% allCohResults               %reshapeCoh
% ------------------------------------------------
% lum 5                 lum 5  | scone 5  |  RG 5
% lum 12                lum 12 | scone 12 |  RG 12
% lum 16                lum 16 | scone 16 |  RG 16
% scone 5
% scone 12
% scone 16
% RG 5
% RG 12
% RG 16

reshapeCoh=mmy_Remove_Harmonic_Powers(reshapeCoh, [inputFreq], maxFreq);
reshapeIncoh=mmy_Remove_Harmonic_Powers(reshapeIncoh, [inputFreq], maxFreq);

% = this is just to make some figures... not important to ANOVA
% make sure the input frequencies are taken out...

figure(106)

reshapeCohFig=reshape(reshapeCoh, [9, 1000, 16]);
reshapeIncohFig=reshape(reshapeIncoh, [9, 1000, 16]);

for i = 1:16
    rsCohFig(:,:,i) = reshapeCohFig(orderForFig,:,i);
    rsIncohFig(:,:,i) = reshapeIncohFig(orderForFig,:,i);
end

for i = 1:length(condTriggers)
    s=subplot(3,3,i);
    bar(rsCohFig(i, myrange));
    h=title(plotTitle{i});
    set(h, 'Visible', 'on');
    linkaxes([s]);
    ylim([0 1])
end

figure(107)

for i = 1:length(condTriggers)
    s=subplot(3,3,i);
    bar(rsIncohFig(i, myrange));
    h=title(plotTitle{i});
    set(h, 'Visible', 'on');
    linkaxes([s]);
    ylim([0 1])
end

% ========================================================================

% gets the RMS of each endo band...
outputCoh=mmy_Mean_Endogenous_Power(reshapeCoh); % 3 colors * 5 avg powers of 5 endo bands * 16 subjects
outputIncoh=mmy_Mean_Endogenous_Power(reshapeIncoh); 

outputCohSPSS=reshape(outputCoh, [12,16])'; % reshape for SPSS ANOVA testing...
outputIncohSPSS=reshape(outputIncoh, [12,16])';

%% Drawing line graphs
%% MY - This bit now computes SEMs...
meanOIncoh=mean(outputIncohSPSS);  
reshapedOutput=reshape(outputIncohSPSS,[16,3,4]);
meanOIncoh=reshape(meanOIncoh,[3 4]);
sdOIncoh=std(outputIncohSPSS);  sdOIncoh=reshape(sdOIncoh,[3 4]);
semIncoh=sdOIncoh/sqrt(size(outputIncohSPSS,1));
% figure(201)
% for endoBand=1:4
%     errorbar(1:3, meanOIncoh(:,endoBand), semIncoh(:,endoBand)); hold on;
% end
% xlim([0.5 3.5]);
% ylabel({'RMS of EEG Power'});
% legend('Alpha: 8-12 Hz', 'Beta: 12-30 Hz', 'Gamma: 30-80 Hz', 'Theta: 4-8 Hz');

figure(200)
hold off
cList=[.1 .1 .1;.7 .3 .3;.1 .5 .7];
for thisColor=1:3
    % NB we reorder the fBands to plot them as monotonically increasing 
    eb(thisColor)=errorbar(1:4, meanOIncoh(thisColor,[4 1 2 3])', semIncoh(thisColor,[4 1 2 3])'); hold on;
    set(eb(thisColor),'LineWidth',2);
    set(eb(thisColor),'Color',cList(thisColor,:));
end
ylabel({'RMS of EEG Power'});
xlim([0 5]);
xlabel('Frequency band');
legend({'Lum','L-M','S'});
set(gca,'XTick',[1 2 3 4 ] );
set(gca,'XTickLabels',{'Theta','Alpha','Beta','Gamma'});
%%
figure(203)
clear g
clear ad
%reshapedOutput=reshape(outputIncoh,[16,3,4]);
ad.dat=reshapedOutput(:,:,[4 1 2 3]);
ad.dat=ad.dat(:);
ad.color=kron(ones(16,1),repmat([1,2,3],1,4));
ad.color=ad.color(:);
ad.band=kron(ones(16*3,1),repmat([1 2 3 4],1,1));
ad.band=ad.band(:);
bandList={'Theta','Alpha','Beta','Gamma',};
for thisName=1:length(ad.band)
    ad.bandName{thisName}=bandList{ad.band(thisName)};
end


g=gramm('x',ad.bandName,'y',ad.dat,'color',ad.color); % This is wrong right now..
g.set_order_options('x',0);
%g.stat_boxplot('notch',1);

g.stat_summary('type','sem','geom','errorbar','dodge',.3)
g.set_stat_options('alpha',.05);
g.set_point_options('base_size',4);
g.set_line_options('base_size',4);
g.set_layout_options('legend_position',[.7 .5 .5 .5]);
g.set_color_options('map',[.4 .4 .4;.8 .5 .5;.2 .6 .7]);
g.draw();
g.export('file_name','powerRMS.pdf','file_type','pdf');
g.export('file_name','powerRMS.svg','file_type','svg');

%% One-Way Repeated Measures ANOVA

% Here, do a repeated measures ANOVA on the allCohResults data at a
% particular frequency.
% We might, later, do a freq x color interaction
% For now we pick 12Hz as a nice 'middle' frequency.

oneF1Freq=[5 12 16];
twoF1Freq=[oneF1Freq*2];

inputFreq=[{'5 Hz', '12 Hz', '16 Hz'}];
compType=[{'Lum v S-cone', 'Lum v L-M', 'S-cone v L-M'}];

%% INCOHERENT RESULTS

thisMainANOVAIncoh=cell(3,99);

for fOfInterest=1:99 % we wont' be interested in anything beyond
    
    % You're essentially running 99 tiny little ANOVAs on all the
    % frequencies.
    
    %fOfInterest=[]; % Hz
    for conditionNo = 1:3 % i.e. the L+M+S, S-iso, and L-M of a single input frequency
        % condition refers to INPUT FREQUENCY
        
        dataToAnalyze=squeeze(allIncohResults([conditionNo, conditionNo+3, conditionNo+6],fOfInterest+1,:))';
        
        t = table(dataToAnalyze(:,1),dataToAnalyze(:,2),dataToAnalyze(:,3),...
            'VariableNames',{'Lum','S','RG'});
        
        colorType=[1 2 3]';
        
        rm = fitrm(t,'Lum-RG ~ 1','WithinDesign', colorType);
        
        ranovatbl = ranova(rm);
        
        mauchlyP = mauchly(rm);
        
        %ranovaary = table2array(ranovatbl);
        
        thisMainANOVAIncoh{conditionNo,fOfInterest,:,:}=ranovatbl;
        
        if mauchlyP.pValue < 0.50 % If it violates the rules of sphericity...
            
            thisPVal(conditionNo,fOfInterest)=ranovatbl.pValueLB(1);
            %disp('Mauchly''s test indicates that assumption of sphericity was violated (p < .05)');
        else
            thisPVal(conditionNo,fOfInterest)=ranovatbl.pValue(1);
        end
        
    end
end

i=1;

for conditionNo=1:3
    for fOfInterest=1:99
        if (thisPVal(conditionNo,fOfInterest) < 0.05) % if it's a sig effect
            thisAnovaStatIncoh(i,1)=thisMainANOVAIncoh{conditionNo,fOfInterest}.DF(1);
            thisAnovaStatIncoh(i,2)=thisMainANOVAIncoh{conditionNo,fOfInterest}.DF(2);
            thisAnovaStatIncoh(i,3)=thisMainANOVAIncoh{conditionNo,fOfInterest}.F(1);
            thisAnovaStatIncoh(i,4)=thisPVal(conditionNo,fOfInterest);
            
            i=i+1;
        end
    end
end


% POST HOC TEST: MULTIPLE COMPARISONS
bonfCount=1;
compCount=1;

compSetA=[1 1 2];
compSetB=[2 3 3];

for fOfInterest=1:99
    
    for conditionNo=1:3
        
        if thisPVal(conditionNo, fOfInterest) < 0.05 % if the effect is significant...
            % Trying the Matlab way for Post Hoc
            
            dataToAnalyze=squeeze(allIncohResults([conditionNo, conditionNo+3, conditionNo+6],fOfInterest+1,:))';
            % you should reshape it in a way that goes V1 V1 V1 V2 V2 V2.. VN.
            dataToAnalyzeVector=reshape(dataToAnalyze, [1 numel(dataToAnalyze)]);
            colorLabel={'Lum','Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum',...
                'Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone',...
                'Scone', 'Scone','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG'};
            
            [~, ~, stats] = anova1(dataToAnalyzeVector, colorLabel, 'off');
            % This gives you a table that matlab understands for running
            % multiple comparisons/post hoc.
            
            bonf=multcompare(stats, [], 'off');
            thisBonfIncoh(:,:,bonfCount)=bonf;
            
            % Level, Level, Lower 95% Confidence, Mean Diff, Upper 95%,
            % P-value
            
            bonfCount=bonfCount+1;
            
            % =========== ACTUAL T TEST ===========
            
            for compNo=1:3
                
                [~,p,~,tStats] = ttest2(dataToAnalyze(:,compSetA(compNo)), dataToAnalyze(:,compSetB(compNo)));
                
                thisTtestIncoh(compNo,1,compCount)=p;
                thisTtestIncoh(compNo,2,compCount)=tStats.tstat;
                thisTtestIncoh(compNo,3,compCount)=tStats.df;
                thisTtestIncoh(compNo,4,compCount)=tStats.sd;
                
            end
            
            compCount=compCount+1;
            
        end
        
    end
    
end

displayFig=0;

if displayFig
    
    figure(1)
    
    for freqVal = 1:3
        title('Incoherent Averages');
        s(freqVal)=subplot(3,1,freqVal);
        bar(-log10(thisPVal(freqVal,:)), 'FaceColor', [0 0.5 0]); hold on;
        line([0 100], [-log10(0.05) -log10(0.05)]); % line that indicates p=0.05
        line([0 100], [-log10(0.01) -log10(0.01)]); % line that indicates p=0.01
    end
    
    hold off;
    linkaxes([s]);
    ylim([0 4])
    
end

saveas(gcf, 'Fig_1_Incoh_100Hz.pdf');

%% COHERENT RESULTS

thisMainANOVACoh=cell(3,99);

for fOfInterest=1:99 % we wont' be interested in anything beyond
    
    % You're essentially running 99 tiny little ANOVAs on all the
    % frequencies.
    
    %fOfInterest=[]; % Hz
    for conditionNo = 1:3 % i.e. the L+M+S, S-iso, and L-M of a single input frequency
        % condition refers to INPUT FREQUENCY
        
        dataToAnalyze=squeeze(allCohResults([conditionNo, conditionNo+3, conditionNo+6],fOfInterest+1,:))';
        
        t = table(dataToAnalyze(:,1),dataToAnalyze(:,2),dataToAnalyze(:,3),...
            'VariableNames',{'Lum','S','RG'});
        
        colorType=[1 2 3]';
        
        rm = fitrm(t,'Lum-RG ~ 1','WithinDesign', colorType);
        
        ranovatbl = ranova(rm);
        
        mauchlyP = mauchly(rm);
        
        thisMainANOVACoh{conditionNo,fOfInterest}=ranovatbl;
        
        if mauchlyP.pValue < 0.50 % If it violates the rules of sphericity...
            thisPVal(conditionNo,fOfInterest)=ranovatbl.pValueLB(1);
            %disp('Mauchly''s test indicates that assumption of sphericity was violated (p < .05)');
            
        else
            thisPVal(conditionNo,fOfInterest)=ranovatbl.pValue(1);
        end
        
    end
end


i=1;

for conditionNo=1:3
    for fOfInterest=1:99
        if (thisPVal(conditionNo,fOfInterest) < 0.05) % if it's a sig effect
            thisAnovaStatCoh(i,1)=thisMainANOVACoh{conditionNo,fOfInterest}.DF(1);
            thisAnovaStatCoh(i,2)=thisMainANOVACoh{conditionNo,fOfInterest}.DF(2);
            thisAnovaStatCoh(i,3)=thisMainANOVACoh{conditionNo,fOfInterest}.F(1);
            thisAnovaStatCoh(i,4)=thisPVal(conditionNo,fOfInterest);
            
            i=i+1;
        end
    end
end

% POST HOC TEST: MULTIPLE COMPARISONS
bonfCount=1;
compCount=1;

compSetA=[1 1 2];
compSetB=[2 3 3];

for fOfInterest=1:99
    
    for conditionNo=1:3
        
        if thisPVal(conditionNo, fOfInterest) < 0.05 % if the effect is significant...
            % Trying the Matlab way for Post Hoc
            dataToAnalyze=squeeze(allCohResults([conditionNo, conditionNo+3, conditionNo+6],fOfInterest+1,:))';
            
            % =========== COMPARISONS ===========
            
            % you should reshape it in a way that goes V1 V1 V1 V2 V2 V2.. VN.
            dataToAnalyzeVector=reshape(dataToAnalyze, [1 numel(dataToAnalyze)]);
            colorLabel={'Lum','Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum','Lum', 'Lum',...
                'Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone','Scone', 'Scone',...
                'Scone', 'Scone','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG','RG', 'RG'};
            
            [~, ~, stats] = anova1(dataToAnalyzeVector, colorLabel, 'off');
            % This gives you a table that matlab understands for running
            % multiple comparisons/post hoc.
            bonf=multcompare(stats, [], 'off');
            thisBonfCoh(:,:,bonfCount)=bonf;
            
            % Level, Level, Lower 95% Confidence, Mean Diff, Upper 95%,
            % P-value
            
            bonfCount=bonfCount+1;
            
            % =========== ACTUAL T TEST ===========
            
            for compNo=1:3
                
                [~,p,~,tStats] = ttest2(dataToAnalyze(:,compSetA(compNo)), dataToAnalyze(:,compSetB(compNo)));
                
                thisTtestCoh(compNo,1,compCount)=p;
                thisTtestCoh(compNo,2,compCount)=tStats.tstat;
                thisTtestCoh(compNo,3,compCount)=tStats.df;
                thisTtestCoh(compNo,4,compCount)=tStats.sd;
                
            end
            
            compCount=compCount+1;
            
        end
        
    end
    
    % Doing a pairwise T test instead..
end





% Bar plot

if displayFig
    
    figure(2)
    
    for freqVal = 1:3
        title('Coherent Averages');
        s(freqVal)=subplot(3,1,freqVal);
        bar(-log10(thisPVal(freqVal,:)), 'FaceColor', [0 0.5 0]); hold on;
        line([0 100], [-log10(0.05) -log10(0.05)]); % line that indicates p=0.05
        line([0 100], [-log10(0.01) -log10(0.01)]); % line that indicates p=0.01
    end
    
    hold off;
    linkaxes([s]);
    ylim([0 4])
    
end

saveas(gcf, 'Fig_2_Coh_100Hz.pdf');

%% Scatterplot??

if displayFig
    
    figure(4)
    
    subplot(3,1,1)
    s1=scatter(1:length(thisPVal),-log10(thisPVal(1,:)), 40, 'MarkerEdgeColor', [0 .5 .5],...
        'MarkerFaceColor', [0 .7 .7], ...
        'LineWidth', 1.5); % GREEN: 5 Hz
    line([0 100], [-log10(0.05) -log10(0.05)]);
    
    subplot(3,1,2)
    s2=scatter(1:length(thisPVal),-log10(thisPVal(2,:)), 40, 'MarkerEdgeColor', [0 0 1],...
        'MarkerFaceColor', [0 0 1], ...
        'LineWidth', 1.5);  % BLUE: 12 HZ
    line([0 100], [-log10(0.05) -log10(0.05)]);
    
    subplot(3,1,3)
    s3=scatter(1:length(thisPVal),-log10(thisPVal(3,:)), 40, 'MarkerEdgeColor', [1 0 0],...
        'MarkerFaceColor', [1 0 0], ...
        'LineWidth', 1.5);  % RED: 16 HZ
    line([0 100], [-log10(0.05) -log10(0.05)]);
    
end

%% Two-way Repeated Measures ANOVA

% figure(1);
% scatter(thisPVal);
%
% for thisCondNo = 1:3
%     for thisFOfInterest = 1:98
%         if thisPVal(thisCondNo,thisFOfInterest) < 0.05
%             scatter(thisCondNo,thisFOfInterest);
%         end
%     end
% end

% bar(-log10(thisPVal))
%
%
% %%
% colorType=[1 2 3]';
% t = table(colorType,dataToAnalyze(:,1),dataToAnalyze(:,2),dataToAnalyze(:,3),'VariableNames',{'colorType','lum', 's', 'rg'});
% meas = table([1 2 3]','VariableNames',{'Measurements'});