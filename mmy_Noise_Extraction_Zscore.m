function zscoreIndices = mmy_Noise_Extraction_Zscore (dataVector, zscoreLimit)

% A function that takes in a vector V, and calculates the difference between
% V(i) and V(i+1), and finds out the zscore. 

% A good way to use this is to input the data vector of the blink
% electrode. You can pick out the indices (blinkIndices) of when blink noise is more than
% one standard deviation away from the mean, indicating a 'big enough'
% blink noise. You can choose to remove these or change them to the mean
% value. 

% Z-score less than 0 represents an element less than the mean. 
% Z-score more than 0 represents an element greater than the mean. 
% 68% of elements have Z score from between -1 and 1;
% 95% of elements have Z score from between -2 and 2;
% 99% of elements have Z score from between -3 and 3.

diffVector = diff(dataVector);
zscoreData = zscore(diffVector);

% zscoreIndices are indices of the data where z score is beyond the denoted
% range.
zscoreIndices = find(abs(zscoreData)>zscoreLimit);



figure(500)
subplot(2,1,1)
plot(dataVector)
title('Data Vector');
subplot(2,1,2)
plot(zscoreData); hold on;
plot(zscoreIndices+1, 0, 'ro'); 
title('Diff vector');
hold off;
end
