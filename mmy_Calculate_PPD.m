function [PPD, viewingAngleDeg] = mmy_Calculate_PPD(viewingDistance, verticalScreenDimensionCm, verticalScreenResolutionPix)
% function to calculate pixels per degree based on given parameters
%all in cm
%written by MY, 11-Sep-2017.

viewingAngleDeg=round(atand((verticalScreenDimensionCm*0.5)/viewingDistance));
PPD=round((verticalScreenResolutionPix*0.5)/viewingAngleDeg);

return
