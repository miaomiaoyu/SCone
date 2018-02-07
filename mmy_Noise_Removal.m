function noiseFilteredData = mmy_Noise_Removal(ftData, range)

meanData=nanmean(abs(ftData),2); 
logData=-log(meanData(range));

g=polyfit(range, logData(:)', 1);
noiseRemoved = logData-range' .* g(1)-g(2);
noiseFilteredData = exp(noiseRemoved);

end