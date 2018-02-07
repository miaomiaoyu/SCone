function params = displayParamsVP
% For ViewPixx (EEG) in B114, University of York, Psychology
% Driven by an Apple Mac Pro
% Ocean Optics Jaz Photospectrometer with fibre optic 
% Calibrated 14-Aug-2017, MY

% Critical parameters
params.numPixels = [1920 1080];
params.dimensions = [52.2 29.3]; %cm
params.distance = 57;
params.frameRate = 120;
params.cmapDepth = 10; % how many levels represented by each gun, used to be 256/8bits; for viewpixx 3D it's 10bit (i.e. 2^10)
params.screenNumber = 0;

% Descriptive parameters
params.computerName = 'tyrion';
params.monitor = 'ViewPixx/3D';
params.card = 'AMD FirePro D500'; % find in systems
params.position = 'PS/B/114';

% parameters which can make programming the display easier
params.flipLR = 0;
params.flipUD = 0;

return
