close all;
clear all

aStart=now;
EEGpath = '/wadelab_shared/Projects/NeuralOscillations/Koniocellular/KCIM_EEG/';
sList={'AN', 'ARW',  'BJ'  ,'CJ' , 'EH' , 'FS' , 'GV' , 'JT' , 'LH',  'MS',  'MTS',  'MY',  'PY',  'RM' , 'SP',  'WD'};

 % We leave one out because they don't have enough trials.

for thisSub=1:length(sList)
    fprintf('\nRunning subject %s : %d of %d\n',sList{thisSub},thisSub,length(sList));
    
    dataOut(thisSub)=arw_analyseERPdataMMLiblinearCleanComplex(sList{thisSub},1,EEGpath);
end
aEnd=now;
fprintf('\nDone\nStarted %s, ended %s\n',datestr(aStart),datestr(aEnd));


for t=1:length(dataOut)
    mLine(t,:,:)=dataOut(t).mKL;
end

%%
gm=squeeze(mean(mLine));
se=squeeze(std(mLine))/sqrt(length(dataOut));


figure(12); hold off;
fc=hsv(15);


xlabel('Time (ms)');
ylabel('Percent correct');
%% --------

figure(12); hold off;
fc=hsv(15);
clear sigPoints;
for t=1:size(mLine,2)
    sigPoints{t}=d_doclustercorr((squeeze(mLine(:,t,:))/100), 1, .5, .05, 1000);
p=shadedErrorBar(1:99,squeeze(gm(t,:)),se(t,:));
    
    set(p.patch,'FaceColor',fc(t,:));
    set(p.mainLine,'Color',fc(t,:))*2;
    set(p.mainLine,'LineWidth',2);
    set(p.patch,'FaceAlpha',.4);
    hold on;
    
    sp=sigPoints{t};
    nSigPoints=length(sp);
    for thisCluster=1:nSigPoints
        thisClusterList=sp{thisCluster};
        for thisSigPointMarker=1:length(thisClusterList)
            thisClustPoint=thisClusterList(thisSigPointMarker);
            sigHand=scatter(thisClustPoint,100-t*2);
            set(sigHand,'MarkerFaceColor',[fc(t,:)]);
            set(sigHand,'MarkerFaceAlpha',.5);
               set(sigHand,'MarkerEdgeAlpha',.0);
        end
    end
    
    
end
xlabel('Time (ms)');
ylabel('Percent correct');

grid on;


% We'd like to add significance to this. One (safeish) way to so this is to
% use the baseline period (0-200ms) as an estimate of the false positive
% rate, then set a cluster threshold based on this.


grid on;
save(datestr(now))
