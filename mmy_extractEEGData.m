function EEG=mmy_extractSubjectData_EEG_ANT(fileName,goodChans,doPlotFlag); % Assumes one rep per subject

%% function EEG=extractSubjectData_EEG_ANT(subjID,repID,dirPath)
% Extracts single subject data block from a particular subject and
% repetition in the given directory path.
% Returns common average as well in case you want it.

% Returns EEG as a struct with the following data:
%
%
%
%
%
%
%
%
%
%
%
%

% Get a list of the valid .cnt files in the directory. Then loop over them



EEG=processcnt3(fileName); % Calls another function called processcnt3.

eventIndex=1;

for thisEvent=1:length(EEG.event)
    q=str2num(EEG.event(thisEvent).type);
    if (~isempty(q))
        EEG.eventList(eventIndex)=q;
        EEG.timeList(eventIndex)=EEG.event(thisEvent).latency;
        eventIndex=eventIndex+1;
    end
    
end

sampleRate=EEG.rate;
EEG.nBinsSec=10;
EEG.dummyBinsSec=1;

% We have now loaded in a single file. We can look thorough the information
% in the trigger to see what the timings are.
% The trigger info is contained in EEG.triggers and the corresponding time
% stamps are in timeList

EEG.timeListZeroed=EEG.timeList-EEG.timeList(1); % Remove the first time from all other values

%In MM's data, condition codes appear as integers like 62, 13 etc. just
%before the first trigger marker ('2')
% The first thing to do is find these condition codes.

%EEG.condCodeIndices=find(EEG.eventList>99 & EEG.eventList<200);