function [output] = mmy_Remove_Harmonic_Powers(meanEEGPowers, inputFreq, maxOutputFreq);

% finds the harmonics according to inputFreq and make the EEG powers 0. 

% input(meanEEGPowers) should be a matrix that is EEG power arranged:
% freq * color * output EEG power (1:N Hz) * no of Subj


nVar1=size(meanEEGPowers,1);
nVar2=size(meanEEGPowers,2);
nOutputHz=size(meanEEGPowers,3);
nSubj=size(meanEEGPowers,4);

for subjNo=1:nSubj
    
    for freqNo=1:length(inputFreq)
        
        harmonicInterval=1:(maxOutputFreq/inputFreq(freqNo)); %harmonics for each input freq
        harmonicFreq=[inputFreq(freqNo)*harmonicInterval,50];

        for i=1:length(harmonicFreq)   
            meanEEGPowers(freqNo,:,harmonicFreq(i)+1,subjNo)=0;
        end
        
    end
    
end

output=meanEEGPowers;