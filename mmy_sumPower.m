function returnedPower=mmy_sumPower(inputSpectrum, freqToSum)
% alphaPower=mmy_sumpower(inputSpectrum, freqToSum)
% Computes the sum of the squared values of amplitudes at the given frequencies.
% e.g. 
% alphaPower=mmy_sumpower(inputSpectrum, [8:12]+1)
% Note the '+1' because Matlab FT starts at 1
% For now this only works for 1D or 2D arrays.
% But - using 'size' and 'reshape' commands it could easily return the
% correct answer for any dimensionality of input...
% 24/03/17, M.Y.
% 07/02/18 ARW - micro edit
%
returnedPower=sum(abs(inputSpectrum([freqToSum] + 1, :)) .^2);

return