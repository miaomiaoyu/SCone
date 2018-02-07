% for making a IM figure

inputFreq = 12;
endoFreq = 8;
noiseLevel=1;
sampleRate = 1000;

inputLinspace = linspace(0, 2*pi*inputFreq, sampleRate);
inputWave=sin(inputLinspace)+rand(1,length(inputLinspace))*noiseLevel;
endoLinspace = linspace(0, 2*pi*endoFreq, sampleRate);
endogenousWave=sin(endoLinspace)+rand(1,length(endoLinspace))*noiseLevel;

endogenousScale=1+endogenousWave;
modulatedWave=inputWave.*endogenousScale;

ftInput=fft(inputWave);
meanInput=abs(ftInput);
ftEndo=fft(endogenousWave);
meanEndo=abs(ftEndo);
ftMod=fft(modulatedWave);
meanMod=abs(ftMod);

figure()
subplot(2,3,1);
plot(inputWave);
xlabel('Amplitude');
ylabel('Sample');

subplot(2,3,2); 
plot(endogenousWave);

subplot(2,3,3);
plot(modulatedWave);

subplot(2,3,4);
bar(meanInput(2:25));
xlim([0 25]);
ylim([0 1000]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(2,3,5);
bar(meanEndo(2:25));
xlim([0 25]);
ylim([0 1000]);

subplot(2,3,6);
bar(meanMod(2:25));
xlim([0 25]);
ylim([0 1000]);



