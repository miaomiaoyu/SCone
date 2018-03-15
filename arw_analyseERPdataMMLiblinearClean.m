function outData=arw_analyseERPdataMMLibLinearClean(varargin)
% function outData=arw_analyseERPdataMMLibLinearClean(varargin)
% script to do basic analysis of EEG data and (optionally) MVPA
% Now requires access to the mex files from the liblinear package (see also libsvm)
% https://www.csie.ntu.edu.tw/~cjlin/liblinear/
% % DHB 11/5/16
% Modified by ARW 05/03/2018
% 07/03/18 : ARW Edited to try and use frequency domain instead of time
% domain.

nbootstrapruns = 1000;
maxFrequency=100; % Look from 1Hz up to this point...
junkBins=1;
goodBins=7;
%
if nargin
    subj = varargin{1};
    processfiles = varargin{2};
    EEGpath=varargin{3};
    
else
    error('No inputs');
end
dofiltering = 1;
getname = 1;
subjname = [];

subjname = char(subjname)

condcodes = [24,59,79,34,83,111,54,131,175]; % Condition codes that you find in the dataset. For MM these are 24,59,79,343,83,111,54, 131, 175
nConds=length(condcodes);
conditioncounter = zeros(1,length(condcodes)); % I think this is blocks x conditions


d = dir(strcat(EEGpath,subj,'/'));
blockcounter = 0;
for n = 1:length(d)
    temp = d(n).name;
    if length(temp)>3
        if temp(end-2:end)=='cnt'
            blockcounter = blockcounter + 1;
            namelist{blockcounter} = temp;
        end
    end
end
fprintf('\nFound %d blocks in %s',blockcounter,strcat(EEGpath,subj));

filename=strcat(EEGpath,subj,'/');
fprintf('\nLoading %s\n',filename);

EEG = processcnt(filename,1);


clear triggertimes triggercond

% For MMY's experiments, codes are 24,59,79,34,83,111,54,131,175 for different
% combinations of colour (Lum, L-M and S) and Frequency (5,12,16 Hz)
%.
%% Get the unique codes
triggercount = 0;
for n = 1:length(EEG.event)
    
    if ~isempty(str2num(EEG.event(n).type))
        thisCode(n)=str2num(EEG.event(n).type);
        if sum(find(condcodes(:)==str2num(EEG.event(n).type)))>0        % now just allow triggers from the list
            triggercount = triggercount + 1;
            triggertimes(triggercount) = EEG.event(n).latency;
            triggercond(triggercount) = str2num(EEG.event(n).type);
        end
    end
end
%%
triggertimes=round(triggertimes); % Round to nearest ms

a = find(ismember(triggercond,condcodes)); % This should compare against condcodes....     5 is a blank. Come to think of it, it would be nice to have this in there as well as a control..
trialtimes = triggertimes(a);
trialconds = triggercond(a);    % remove the first digit so codes are 1-4

starttrial = 1;
ntrials = length(trialtimes);

% PRe-allocate a big array of zeros that will contain the individual bin
% data for each condition type. 
% This will be 99 (maxFrequency-1) x nGoodBins * number of instances of
% that trial type (15 in this case) x conditions
nInstancesOfEachCondition=15;
nSensors=66;
%%
rawFFTData=zeros(nSensors,maxFrequency-1,goodBins,nInstancesOfEachCondition,nConds);
tic
for trial = starttrial:ntrials
    
    currenttrial = trialconds(trial);
    trialIndex=find(condcodes==currenttrial);
    
    conditioncounter(trialIndex) = conditioncounter(trialIndex) + 1;
     
  
    temp = EEG.data(:,(trialtimes(trial)+(junkBins*EEG.rate)+1):(trialtimes(trial)+((junkBins+goodBins)*EEG.rate)));
    % Pull out all the good data bins
    % Reshape them:
    tempR=reshape(temp,66,EEG.rate,goodBins);
    % Now compute the ft
    fTempR=fft(tempR,[],2);
    %Chop at the right frequency
    fTempR_Chopped=fTempR(:,2:maxFrequency,:); % For now this is still complex....
    
    
    % Thhis little bit here.... you need to add the reshaped matrix you
    % just made into a big array (alltrials) in the correct place. This
    % will be 150 * 9 * 66 * 100 :  bins x conds x channels x frequency points 
    %alltrials(currenttrial,conditioncounter(currenttrial),:,:) = resampData'; 
    rawFFTData(:,:,:,conditioncounter(trialIndex),trialIndex)=fTempR_Chopped;
    
    
     
end
toc
%%
%

rawFFTData=reshape(rawFFTData,[nSensors,maxFrequency-1,goodBins*nInstancesOfEachCondition,nConds]);

% 
% for n = 1:EEG.nchan
%     if EEG.chanlocs(n).labels(1:2)=='Oz'
%         targetchannelnumber(1) = n;
%     end
%     if EEG.chanlocs(n).labels(1:2)=='O1'
%         targetchannelnumber(2) = n;
%     end
%     if EEG.chanlocs(n).labels(1:2)=='O2'
%         targetchannelnumber(3) = n;
%     end
%     
%     if length(EEG.chanlocs(n).labels)>2
%         if EEG.chanlocs(n).labels(1:3)=='POz'
%             targetchannelnumber(4) = n;
%         end
%         if EEG.chanlocs(n).labels(1:3)=='PO4'
%             targetchannelnumber(5) = n;
%         end
%         if EEG.chanlocs(n).labels(1:3)=='PO6'
%             targetchannelnumber(6) = n;
%         end
%         if EEG.chanlocs(n).labels(1:3)=='PO8'
%             targetchannelnumber(7) = n;
%         end
%         if EEG.chanlocs(n).labels(1:3)=='PO3'
%             targetchannelnumber(8) = n;
%         end
%         if EEG.chanlocs(n).labels(1:3)=='PO5'
%             targetchannelnumber(9) = n;
%         end
%         if EEG.chanlocs(n).labels(1:3)=='PO7'
%             targetchannelnumber(10) = n;
%         end
%     end
% end




% % subtract out pre-trial baseline for all trials
% fprintf('\nThis is the size of ''s''\n');
% 
% s = size(alltrials)
% 
% for cond = 1:s(1)
%     for trial = 1:s(2)
%         for ch = 1:s(3)
%             temp = squeeze(alltrials(cond,trial,ch,:));
%             temp = temp - mean(temp(1:50));
%             alltrials(cond,trial,ch,:) = temp;
%         end
%     end
%     
% end
% %%
% figure(11);
% meanTS=squeeze(mean(alltrials,2));
% semTS=squeeze(std(alltrials,[],2))/sqrt(220);
% errorbar(squeeze(meanTS(:,32,:))',squeeze(semTS(:,32,:))');

%%
tic
nsamplespermean = 5;       % must divide into 105 as an integer

%%% Stoppinghere 15/3/2018. Todo: Convert fft stuff to abs (or phase). 
%% Probably reformat that big FFT array into the same shape as alltrials to make it fit below: (conds x trials x channels x freq)
% 1. CD towards
% 2. CD away
% 3. IOVD towards
% 4. IOVD away
complistA = [1 3 1 2 1]; % Comparisons. You compare one thing from a with one thing from b.So first comp is 1 v 2, then 3 v 4 etc...

complistB = [2 4 3 4 5];
allmvpa = zeros(length(complistA),nResampPoints);  % matrix to store the MVPA results



alltrialsN=alltrials(:,:,1:64,:); % I don't know. Eye channels?



clear datameanvectA;
clear datameanvectB;
clear allScorePred;
clear allKFoldLoss;
clear allMeanPred;
clear allStdPred;
nComparisons=length(compListA);

allKFoldLoss=zeros(nComparisons,nbootstrapruns,nResampPoints);



for comp = 1:9      % three comparisons
    comp % Display the current condition
    totalsamples(1) = 210;           % however many trials we have per condition. Because we want to block these into groups, we round to some easily-factored number
    totalsamples(2) = 210;
    M = totalsamples/nsamplespermean; % This is however many individual points we try to classify from this class each time. Each point is generate by computing a mean of dozens of original time series.
    labels=[ones(M(1),1);ones(M(2),1)*-1];
    tic % Time each condition...
    
    
    for runno = 1:nbootstrapruns % Repeat the sampling over a large number of bootstrapped resamples of different averaged sets
        Aindices = randperm(totalsamples(1)); % Randomly permute the set of sample indices.
        Bindices = randperm(totalsamples(2));
        AdataCond = squeeze(alltrialsN(complistA(comp),:,:,:)); % Then pick out a set relevant to the conditions we are looking at right now.
        BdataCond = squeeze(alltrialsN(complistB(comp),:,:,:));
        
        datameanvectACond=(squeeze(mean(reshape(AdataCond(Aindices,:),nsamplespermean,M(1),64,nResampPoints))));
        datameanvectBCond=(squeeze(mean(reshape(BdataCond(Aindices,:),nsamplespermean,M(2),64,nResampPoints))));
        %c=sprintf('a=train(l,d,''-c 1 -v 5 -q'')'); % This sets up the call to 'train' (which is part of the nice mex file liblinear). Note we make it a string so that we can call it with 'eval' so that we can catch output to console. I might edit the mex file to make this unnecessary..
        % Note also - you might feel tempted to use the -n option to
        % increase the number of threads on a multi-cpu machine. As far
        % as I can tell, for these types of small instance, high
        % dimensionality datasets, this actually slows things down. It's
        % very fast anyway....
        
        l=double(labels); % Everything going into liblinear must be a double
        da=double(datameanvectACond(:,:,:));
        db=double(datameanvectBCond(:,:,:));
        for t = 1:nResampPoints % Time points
            
            
            
            d=sparse(([da(:,:,t);db(:,:,t)])); % This is the way liblinear likes its inputs
            
            a=train(l,d,'-c 1 -v 5 -q');
            %  q=evalc(c); % Doing this to avoid message to console.
            %This calls liblinear's 'train' routine. Because we've asked to do kfold validation, this will simply return a single number in 'a' which is the classifier accuracy. Random chance is 50,
            
            allKFoldLoss(comp,runno,t)=a; % Keep track of each bootstrap iteration's accuracy
            
        end % next t
        
    end % next boot
    toc
    
end % Next cond

% All we really need to return is the mean over bootstraps and the std.
mKL=squeeze(mean(allKFoldLoss,2));
sKL=squeeze(std(allKFoldLoss,[],2));

% Some plotting. Really (for a function) we don't want this.
akf=shiftdim(allKFoldLoss,1);



figure(10);hold off;
fc=[.3 0 0;0 .3 0;0 0 .3; .1 .1 .1];

for t=1:size(akf,3)
    
    p=shadedErrorBar([],mean(akf(:,:,t)),std(akf(:,:,3)));
    set(p.patch,'FaceColor',fc(t,:));
    set(p.mainLine,'Color',fc(t,:))*2;
    set(p.mainLine,'LineWidth',2);
    set(p.patch,'FaceAlpha',.4);
    hold on
end


outData.mKL=mKL;
outData.sKL=sKL;
outData.nbootstrapruns=nbootstrapruns;
outData.subj=subj;

return
%%
end








%--------------------------------------------------------------------------
function EEG = processcnt(dirpath,block)

d = dir(dirpath);
counter = 0;
for n = 1:length(d)
    temp = d(n).name;
    if length(temp)>3
        if temp(end-2:end)=='cnt'
            counter = counter + 1;
            namelist{counter} = temp;
        end
    end
end

[chanlocs, lay] = d_loadmontage;

filename = strcat(dirpath,namelist{block});

tempeeg = read_eep_cnt(filename,1,2);

EEG = read_eep_cnt(filename,1,tempeeg.nsample);
EEG.data = single(EEG.data);        % halves the disc space required - this is what EEG lab does automatically and doesn't make any difference as far as I can tell

filename = strcat(filename(1:end-3),'trg');

tempevent = read_eep_trg(filename);

for ev = 1:length(tempevent)
    EEG.event(ev).type = tempevent(ev).code;
    EEG.event(ev).latency = tempevent(ev).time;
end
EEG.chanlocs = chanlocs;
EEG.lay = lay;


end
%--------------------------------------------------------------------------
function channelmappings = getchannelmappings(hmod,eegchans)

for n = 1:66
    labeltofind = upper(eegchans(n).labels);
    for m = 1:66
        labeltocheck = hmod.label{m};
        minlength = min(length(labeltofind), length(labeltocheck));
        if labeltofind(1:minlength)==labeltocheck(1:minlength)
            channelmappings(n) = m;
        end
    end
end

end
%--------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% ends an array 'mask' that is 1 inside the circular window, shading to zero outside
% W, H are the width and height of the whole (rectangular or square) array, in pixels
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is the smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies this diameter (in relative units, range 0 -> 1)
% MAG, 27.2.04

%soft window parameters
if nargin<3, D = 0.9; end % sets default diameter to 0.9
radius = min(W*D/2,H*D/2);% radius in pixels
blur = 2*(min(W/2,H/2) - radius);  % blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
%image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% make circular soft window
mask = single((xx.*xx + yy.*yy) < radius^2); % logical 0 outside the circle,1 inside it
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask/max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%--------------------------------------------------------------------------------------------------