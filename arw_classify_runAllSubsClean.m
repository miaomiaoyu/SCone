
aStart=now;
EEGpath = '/wadelab_shared/data/MiD_EEG_Classify/';

sList={'AW','BP','DL','FS','GM','IJ','IS','JS','MW','RM'}; % We leave one out because they don't have enough trials.

for thisSub=1:length(sList)
    fprintf('\nRunning subject %s : %d of %d\n',sList{thisSub},thisSub,length(sList));
    
    dataOut(thisSub)=analyseERPdataFedLiblinearClean(sList{thisSub},1,EEGpath);
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
fc=[.3 0 0;0 .3 0;0 0 .3; .1 .1 .1];

for t=1:4
p=shadedErrorBar(linspace(1,1000,125),squeeze(gm(t,:)),se(t,:));
  set(p.patch,'FaceColor',fc(t,:));
    set(p.mainLine,'Color',fc(t,:))*2;
    set(p.mainLine,'LineWidth',2);
    set(p.patch,'FaceAlpha',.4);
hold on;

end
xlabel('Time (ms)');
ylabel('Percent correct');

grid on;
