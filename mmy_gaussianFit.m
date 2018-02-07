function [data_to_fit,fit_coefs,gof, recovered_series,fft_data,log_data,noise_removed] = mmy_gaussianFit(data,linenoise,rangeOI)
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

if isempty(linenoise)
    linenoise=50.0;
    fprintf('line noise taken to be 50 Hz');
end

fft_data=fft(data); % calculate the fft of binned data set
    
mean_data=nanmean(abs(fft_data),2); 
% compute incoherent mean of this fft; 2 means 2nd dimension - per sec.
                                    
log_data=-log(mean_data(rangeOI)); % we're only interested in 2:200 Hz really
g=polyfit(rangeOI,log_data(:)',1); % g(1) gives us the slope, g(2) gives the intercept.

noise_removed=log_data-(rangeOI)'.*g(1)-g(2); % 1/f noise
recovered_series=exp(-noise_removed); % data without the 1/f noise
data_to_fit=recovered_series(rangeOI); % pick out range of interest

nm=nanmean(data_to_fit);

data_to_fit(linenoise./2+1)=nm;data_to_fit(linenoise+1)=nm;%linenoise removal
data_to_fit(find(data_to_fit<0))=nm;

%what is this for
%data_to_fit=data_to_fit./sum(data_to_fit.^2); % this creates NaN when it's 0./0
%data_to_fit(isnan(data_to_fit))=nm; % remove NaN
data_to_fit=data_to_fit-1;
i=1:length(data_to_fit); % no of Hz we're interested in.
[fit_coefs,gof]=matlab_cftool_1(i(:),data_to_fit(:));
%fit(i(:),data_to_fit(:),'gauss4'); % gauss4 copes w if there are multiple peaks

% coeffvals=coeffvalues(fitted_curve)';
% a1=coeffvals(1);b1=coeffvals(2);c1=coeffvals(3);
% a2=coeffvals(4);b2=coeffvals(5);c2=coeffvals(6);

% a - 
% b - x coord of peaks (freq)

% figure;
% subplot(3,1,1);
% plot(log_data);
% title('FFT data');
% subplot(3,1,2);
% plot(-noise_removed);
% title('1/f noise');
% subplot(3,1,3);
% plot(recovered_series);
% title('Data after removing 1/f noise');


% subplot(4,1,4);
% plot(abs(fft_data(2:200)));
% title('abs fft data');
