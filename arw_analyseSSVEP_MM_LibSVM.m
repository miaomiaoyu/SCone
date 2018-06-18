function outData=arw_analyseSSVEP_MM_LibSVM(varargin)
% function outData=arw_analyseERPdataMMLibLinearClean(varargin)
%   subj = varargin{1};
%    processfiles = varargin{2};
%    EEGpath=varargin{3};
%    complexFlag=varargin{4};
%    nbootstrapruns=varargin{5};
%     compArray=varargin{6};
% script to do basic analysis of EEG data and (optionally) MVPA
% Now requires access to the mex files from the liblinear package (see also libsvm)
% https://www.csie.ntu.edu.tw/~cjlin/liblinear/
% % DHB 11/5/16
% Modified by ARW 05/03/2018
% 07/03/18 : ARW Edited to try and use frequency domain instead of time
% domain.

maxFrequency=150; % Look from 1Hz up to this point...
junkBins=1;
goodBins=11;
% c
if nargin
    subj = varargin{1};
    processfiles = varargin{2};
    EEGpath=varargin{3};
    complexFlag=varargin{4};
    nbootstrapruns=varargin{5};
    compArray=varargin{6};
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
rawFFTData=zeros(nConds,nInstancesOfEachCondition,goodBins,nSensors,maxFrequency-1);%,nInstancesOfEachCondition,nConds);
%cond x trial x channel x freq
tic
for trial = starttrial:ntrials
    
    currenttrial = trialconds(trial);
    trialIndex=find(condcodes==currenttrial);
    
    conditioncounter(trialIndex) = conditioncounter(trialIndex) + 1;
    
    temp = EEG.data(:,(trialtimes(trial)+(junkBins*EEG.rate)+1):(trialtimes(trial)+((junkBins+goodBins)*EEG.rate)));
    % Pull out all the good data bins
    % Reshape them:
    tempR=reshape(temp,[66,EEG.rate,goodBins]);
    % Now compute the ft
    fTempR=fft(tempR,[],2);
    %Chop at the right frequency
    fTempR_Chopped=fTempR(:,2:maxFrequency,:); % For now this is still complex....
    
    
    % Thhis little bit here.... you need to add the reshaped matrix you
    % just made into a big array (alltrials) in the correct place. This
    % will be 150 * 9 * 66 * 100 :  bins x conds x channels x frequency points 
    %alltrials(currenttrial,conditioncounter(currenttrial),:,:) = resampData'; 
    %rawFFTData(:,:,:,conditioncounter(trialIndex),trialIndex)=fTempR_Chopped;
    rawFFTData(trialIndex,conditioncounter(trialIndex),:,:,:)=shiftdim(fTempR_Chopped,2);
     

end
toc
%%

rawFFTData=reshape(rawFFTData,[nConds,goodBins*nInstancesOfEachCondition,nSensors,maxFrequency-1]);


%%
tic
nsamplespermean = 5;       % must divide into 105 as an integer. We will chop the fft  into lumps of 5 bins, average them and compute classificaiotn.

% MMY- Here let us know what the condition codes mean...

complistA = compArray(1,:);

complistB = compArray(2,:);



nComparisons=length(complistA);
allmvpa = zeros(nComparisons,maxFrequency-1);  % matrix to store the MVPA results for each bootstrap

alltrialsN=(rawFFTData(:,:,1:64,:)); % Remove eye channels and, for now, take just the abs of each complex number. Later we will use both real and imag parts by assiging twice as many channels and using half for real and half for imag



clear datameanvectA;
clear datameanvectB;
clear allScorePred;
clear allKFoldLoss;
clear allMeanPred;
clear allStdPred;

allKFoldLoss=zeros(nComparisons,nbootstrapruns,maxFrequency-1);

nBinsTotal=size(alltrialsN,2);

for comp = 1:nComparisons      % three comparisons
    comp % Display the current condition
    totalsamples= [nBinsTotal nBinsTotal];           % however many trials we have per condition. Because we want to block these into groups, we round to some easily-factored number
    M = totalsamples/nsamplespermean; % This is however many individual points we try to classify from this class each time. Each point is generate by computing a mean of dozens of original time series.
    labels=[ones(M(1),1);ones(M(2),1)*-1];
    tic % Time each condition...
    
    
    parfor runno = 1:nbootstrapruns % Repeat the sampling over a large number of bootstrapped resamples of different averaged sets
        Aindices = randperm(totalsamples(1)); % Randomly permute the set of sample indices.
        Bindices = randperm(totalsamples(2));
        AdataCond = squeeze(alltrialsN(complistA(comp),:,:,:)); % Then pick out a set relevant to the conditions we are looking at right now.
        BdataCond = squeeze(alltrialsN(complistB(comp),:,:,:));
        
        datameanvectACond=(squeeze(mean(reshape(AdataCond(Aindices,:),nsamplespermean,M(1),64,(maxFrequency-1)))));
        datameanvectBCond=(squeeze(mean(reshape(BdataCond(Aindices,:),nsamplespermean,M(2),64,(maxFrequency-1)))));
        %c=sprintf('a=train(l,d,''-c 1 -v 5 -q'')'); % This sets up the call to 'train' (which is part of the nice mex file liblinear). Note we make it a string so that we can call it with 'eval' so that we can catch output to console. I might edit the mex file to make this unnecessary..
        % Note also - you might feel tempted to use the -n option to
        % increase the number of threads on a multi-cpu machine. As far
        % as I can tell, for these types of small instance, high
        % dimensionality datasets, this actually slows things down. It's
        % very fast anyway....
        
        l=double(labels); % Everything going into liblinear must be a double
        if (complexFlag==1)
            da=double(cat(2,real(datameanvectACond),imag(datameanvectACond)));
            db=double(cat(2,real(datameanvectBCond),imag(datameanvectBCond)));
        else
            da=double(abs(datameanvectACond));
            db=double(abs(datameanvectBCond));
        end
        

        
        % Here we allow complex numbers by essentially doubling the number
        % of electrodes
        
        for t = 1:(maxFrequency-1) % Time points
            
            
            
            d=(([da(:,:,t);db(:,:,t)])); % This is the way liblinear likes its inputs
            d=zscore(d,[],2)
            a=libsvmtrain(l,d,'-c 1 -v 5 -q');
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
fc=hsv(15);


for t=1:size(akf,3)
    
    p=shadedErrorBar([],mean(akf(:,:,t)),std(akf(:,:,t)));
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
