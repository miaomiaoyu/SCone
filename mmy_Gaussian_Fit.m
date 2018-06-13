function [dataToFit, fitCoefs, GOF, recoveredSeries, fftData, logData, noiseRemoved]...
    = mmy_Gaussian_Fit (data, lineNoise, rangeOfInterest)
% mean_data,log_data,recovered_series removed 
% mmy_gaussianfit     takes in (1) data: matrix of nSamplesPerSecond by 
%                     nSeconds (2) lineNoise: in Hz (for noise removal) and 
%                     (3) HzRange (of what you're interested in). Returns 
%                     the peak amplitude, frequency and goodness of fit based
%                     on trying to fit a Gaussian curve to the data. Can be
%                     used for electrophysiology data for rodents and
%                     humans.
%
% Notes on Fitting:
%
% The peak amplitude and frequency for each set of data is obtained throgh
% trying to fit broad-peaked function that looks a bit like a Gaussian. You
% can then read off the properties about that fit, like the centre, the
% height, the width and the goodness of fit. 
% 
% To do this, you have to first do two things: (1) get rid of
% harmonics of 50Hz - (just set 25, 50, 100, 150 to 0). (2) Fit.
%
% With fitting, what you're really trying to do is remove 1/f
% noise. A good way to do this is the log transform the y axis and
% fit a linear function to the data because (log(1/x)) is -x.
%         
%
% see also: mmy_bioAnalyzePower, mmy_computeBinned, importaxo
% written by MY, 25-Jul-2017 ---------------------------------------------

if isempty(lineNoise)
    lineNoise=50.0;
    fprintf('Line Noise assumed to be 50 Hz');
end

fftData=fft(data); % calculate the fft of binned data set
    
meanData=nanmean(abs(fftData),2); 
% compute incoherent mean of this fft; 2 means 2nd dimension - per sec.
                                    
logData=-log(meanData(2:200)); % we're only interested in 2:200 Hz really
g=polyfit(2:200,logData(:)',1); % g(1) gives us the slope, g(2) gives the intercept.

noiseRemoved=logData-(2:200)'.*g(1)-g(2); % 1/f noise
recoveredSeries=exp(-noiseRemoved); % data without the 1/f noise
dataToFit=recoveredSeries(rangeOfInterest); % pick out range of interest

nm=nanmean(dataToFit);

dataToFit(lineNoise./2+1)=nm;
dataToFit(lineNoise+1)=nm; %linenoise removal
dataToFit(find(dataToFit<0))=nm;

%what is this for
%data_to_fit=data_to_fit./sum(data_to_fit.^2); % this creates NaN when it's 0./0
%data_to_fit(isnan(data_to_fit))=nm; % remove NaN
dataToFit=dataToFit-1;
i=1:length(dataToFit); % no of Hz we're interested in.
[fitCoefs, GOF]=mmy_Curve_Fit_Tool(i(:), dataToFit(:));
%fit(i(:),data_to_fit(:),'gauss4'); % gauss4 copes w if there are multiple peaks

% coeffvals=coeffvalues(fitted_curve)';
% a1=coeffvals(1);b1=coeffvals(2);c1=coeffvals(3);
% a2=coeffvals(4);b2=coeffvals(5);c2=coeffvals(6);

% a - 
% b - x coord of peaks (freq)

% figure ()
% subplot(4,1,1);
% plot(logData(rangeOfInterest));
% title('FFT data');
% subplot(4,1,2);
% plot(-noiseRemoved(rangeOfInterest));
% title('1/f noise');
% subplot(4,1,3);
% plot(recoveredSeries(rangeOfInterest));
% title('Data after removing 1/f noise');
% subplot(4,1,4);
% plot(abs(fftData((rangeOfInterest))));
% title('abs fft data');
