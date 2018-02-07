function [timingData]=mmy_extractEvents1(EEG)
% function [timingData]=mmy_extractEvents(EEG)
% Accepts a block of EEG data that comes from mmy_extractSubjectData_EEG_ANT
% and returns the identify and start of the unique events (codes 100-199) and, for each
% event, the start of the actual data triggers (codes '2') == for us it's
% 13? I think...
% Obviously this only works if your EEGdata follow these rules 
% 100 <= condCode <200    WHAT IS THIS  - 
% Trigger = 2

% Build our own event list - 
nEventsTotal=length(EEG.eventList); % This is a typical line of code that can be used to 
                                   % comb through a for loop.
for thisEvent=1:nEventsTotal
    a=str2num(EEG.eventList(thisEvent).code); % Convert the string to a number?
    if(isempty(a))
        
        evCode(thisEvent)=0; % change the evCode to 0 if a is empty, but if 'a' is a number,
                             % then add 'a' as the evCode.
    else
        evCode(thisEvent)=a;
    end
    
end
evCode(isempty(evCode))=0;


% UIDIndices=find(evCode>100 & evCode<200); Ours is 13 and 65?

UIDIndices=find(evCode<100);
timingData.allConds=evCode(UIDIndices);

triggerPointIndices=find(evCode==1); % it used to be evCode==2

% QUESTION: what is the EEG. things - are they functions written by us,
% toolbox, module etc.?

UIDs=EEG.eventList(UIDIndices); 
nEvents=length(UIDs);

for thisEvent=1:nEvents
    
    % Make a note of each unique event and then list the time points of the
    % individual triggers that follow that event (up to the next event or
    % '255'
    
    timingData.event(thisEvent).code=UIDs(thisEvent);
    
    
    if (thisEvent<nEvents) % All but the last condition are bounded by another, later condition
                           % If the current event is not the last event,
                           % attach the current index to it (greater than
                           % the current UIDIndice, but smaller than the next one).
                           
        triggerIndicesThisCond=find(triggerPointIndices>UIDIndices(thisEvent) & triggerPointIndices<(UIDIndices(thisEvent+1)));
    else
        triggerIndicesThisCond=find(triggerPointIndices>UIDIndices(thisEvent)); % Unless it's the last event, in which
                                                                                % case just give it the current index?
    end
    
    timingData.triggers(thisEvent).eventIndex=triggerPointIndices(triggerIndicesThisCond);
    for thisOffset=1:length(triggerIndicesThisCond)
        
     timingData.triggers(thisEvent).offsets(thisOffset)=EEG.triggers(triggerPointIndices(triggerIndicesThisCond(thisOffset))).offset;
    end
    
    timingData.triggers(thisEvent).nEvents=length(triggerIndicesThisCond);
    
    
end

