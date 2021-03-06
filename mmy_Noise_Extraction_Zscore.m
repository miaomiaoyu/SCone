function [trashBits, zsIndices, sumBPS] = mmy_Noise_Extraction_Zscore(dataVector, zsLimit, nDataPoints, binSize, nBins, displayFig)

% A function that takes in a vector V, and calculates the difference between
% V(i) and V(i+1), and finds out the zscore.

% A good way to use this is to input the data vector of the blink
% electrode. You can pick out the indices (blinkIndices) of when blink noise is more than
% one standard deviation away from the mean, indicating a 'big enough'
% blink noise. You can choose to remove these or change them to the mean
% value.

% It then takes the vector zsIndices and calculates how
% 'dense' the cluster of data points with zscore beyond zsLimit is for each
% bin. binSize can be sampling rate, to get a score for each second.
% However, I've noted that a better binSize might be 500ms instead, as an
% average blink lasts 300-400ms (and so I don't want to throw away a whole
% second of data). This means that nBins = totalSeconds * 2 (as you get 1
% sec every 2 * 500 samples). 

% Z-score less than 0 represents an element less than the mean.
% Z-score more than 0 represents an element greater than the mean.
% 68% of elements have Z score from between -1 and 1;
% 95% of elements have Z score from between -2 and 2;
% 99% of elements have Z score from between -3 and 3.

% 3.5 seems to be a good number for zsLimit for blink noises. 

% In all: 
% dataVector -> a vector containing all data from blink electrode
% zsLimit -> 3.5
% nDataPoints -> 

%dataVector = blinkPoints;
diffVector = diff(dataVector);
zscoreData = zscore(diffVector);

% zscoreIndices are indices of the data where z score is beyond the denoted
% range.
zsIndices = find(abs(zscoreData)>zsLimit);
%displayFig = 0; nBins = exptEndTime;

r = round(rand(1) * 100000);

markerLineX = 1:binSize-1:10000;

if displayFig
    
    figure(500)
    
    subplot(2,1,1)
    plot(dataVector(r:r+10000-1)); hold on;
    h=fileExchange_Vline(markerLineX, 'g'); 
    
    xlim([0 10000]);
    title('Data Vector');
    
    hold off;
    
    subplot(2,1,2)
    plot(zscoreData(r:r+10000-1)); hold on;
    plot(zsIndices, 0, 'ro'); hold on; % 15th 
    h;

    xlim([0 10000]);
    title('Diff Vector');
   
    hold off;
    
end

bps = zeros(nDataPoints, 1);
bps(zsIndices) = 1;
bps = reshape(bps, binSize, nBins);
sumBPS = sum(bps);

% Note: sumBPS gives you a vector that denotes how many data points beyond
% the zsLimit there are in each bin. You can choose to throw

trashIndices = find(sumBPS > 10);

for i = 1:length(trashIndices)
    trashPoints = trashIndices(i)*binSize;
    trashBits(i,:) = (trashPoints - binSize + 1) : trashPoints;
end

%figure(500)
%subplot(2,1,1)
%plot(dataVector)
%title('Data Vector');
%subplot(2,1,2)
%plot(zscoreData); hold on;
%plot(zscoreIndices+1, 0, 'ro'); 
%title('Diff vector');
%hold off;
