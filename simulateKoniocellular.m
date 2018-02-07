% Simulate what we might find in the KC IM experiment


inputFreq=12; % In Hz.
endogenousFreqBase=8.7; % In Hz

noiseLevel=2; % 1/f pink noise. Or white noise to begin with

sampleRate=1000; % Hz

nBins=10; % How many bins to sample

couplingLevel=2; % This determines how much the endogenous freq modulates the input. For Luminance it will be zero perhaps. For KC it might be a lot.
sumPowerSpectrum=zeros(1000,1);
sumCohAv=zeros(1000,1);

for thisStep=1:100
 endogenousFreq=endogenousFreqBase+rand(1)*2-1;
 
% 1: Simulate two waveforms. One is the input. One is the endogenous freq.
inputSupport=linspace(0,2*pi*inputFreq*nBins,sampleRate*nBins);%+rand(1)*360;

inputWave=sin(inputSupport)+rand(1,length(inputSupport))*noiseLevel;

endogenousSupport=linspace(0,2*pi*endogenousFreq*nBins,sampleRate*nBins)+rand(1)*360;

endogenousWave=sin(endogenousSupport)+rand(1,length(inputSupport))*noiseLevel;



% Now modulate the input by the endogenous thing. Here we will multiply it
% by some offset versoin of the wave
endogenousScale=1+endogenousWave.*couplingLevel;

modulatedWave=(inputWave.*endogenousScale).^2;
subplot(2,1,2);
plot(modulatedWave); % This modulated wave is basially your data....
                    % you have 10 seconds of 1000 samples. so it's 10000
                    % samples altogether. THis is where you jump in
                    % currently. 

% Now look at the FTs

% In real life we will do this: Chop the waveform into 1s bins, compute the
% FT, average incoherently.
finalWave=modulatedWave;% + modulatedWave.^2;

reshapedWave=reshape(finalWave,sampleRate,nBins);
ftModulateWave=fft(reshapedWave);
meanPowerSpectrum=squeeze(mean(abs(ftModulateWave),2)); % Average across the bins .  Remember to average abs rather than the raw FT
meanCohAv=squeeze(mean(ftModulateWave,2));

sumPowerSpectrum=sumPowerSpectrum+meanPowerSpectrum;
sumCohAv=sumCohAv+meanCohAv;

end
figure(1);
subplot(2,1,1);
hold off;
plot(inputWave,'k');
hold on;
plot(endogenousWave,'r');
legend({'Input','Endogenous'});
figure(2);
bar(meanPowerSpectrum(2:25)); % Ignore the first component (0 power - DC offset). Look at the frequewncies from 1Hz to 100Hz.

figure(20);
subplot(2,1,1);
h1=bar(abs(sumCohAv(2:26)));
title('Coherent average');
subplot(2,1,2);
h2=bar(sumPowerSpectrum(2:26));
title('Incoherent average');
h1Bars=get(h1,'Children');
set(h1Bars(12),'FaceColor',[1 0 0]);

