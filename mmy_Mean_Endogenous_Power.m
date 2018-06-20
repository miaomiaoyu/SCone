function [output]=mmy_Mean_Endogenous_Power(meanEEGPowers)

% input should be a matrix that is EEG power arranged:
% var1 * var2 * output EEG power (1:N Hz) * no of Subj

% e.g. for Chromatic SSVEP:
% meanEEGPowersAllSubj = 3 freq * 3 col * 1000 Hz * 16 subj
% freqRange can be any range you want, but in this case we're looking at
% the endogenous rhythms, so they'll correspond to a frequency band.

nFreq=size(meanEEGPowers,1);
nColor=size(meanEEGPowers,2);
nOutputHz=size(meanEEGPowers,3);
nSubj=size(meanEEGPowers,4);

alphaHz=8:12;
betaHz=12:30;
gammaHz=30:80;
deltaHz=1:4;
thetaHz=4:8;

endoHz={alphaHz,betaHz,gammaHz,deltaHz,thetaHz};

for subjNo=1:nSubj
    for typeNo=1:length(endoHz)
        for colorNo=1:nColor
            meanEEGPowerRange=meanEEGPowers(:, colorNo, endoHz{typeNo}+1, subjNo);
            
            meanEEGPowerRangeColor=squeeze(meanEEGPowerRange); % squeeze out the redundant Freq column.
            
            % meanEEGPowerRangeColor is 3 colors * endo freq band (4 for
            % delta, 5 for alpha)
                        
            meanEEGPowerPerBand=sqrt(mean((abs(meanEEGPowerRangeColor).^2),2)); % mean across each frequency contained within an endogenous frequency band... RMS power
            
            % This should give you an averaged EEG power 
            %meanEEGPowerRange=squeeze(meanEEGPowerRange);
        end
        
        endoPower(:,typeNo)=meanEEGPowerPerBand; % 3 colors * averaged endo freq band .. so 3 color by 5 endo bands
        
    end
    
    output(:,:,subjNo)=endoPower;

end
