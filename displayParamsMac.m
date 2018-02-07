function params = displayParamsMac
% For ViewPixx (EEG) in B114, University of York, Psychology
% Driven by an Apple Mac Pro
% Ocean Optics Jaz Photospectrometer with fibre optic 
% Calibrated 14-Aug-2017, MY

% Critical parameters
params.numPixels = [1440 900];
params.dimensions = [28.5 18.3]; %cm
params.distance = 55;
params.frameRate = 60;
params.cmapDepth = 10; % how many levels represented by each gun, used to be 256/8bits; for viewpixx 3D it's 10bit (i.e. 2^10)
params.screenNumber = 0;

% Descriptive parameters
params.computerName = 'Eleanor Rose';
params.monitor = 'Macbook Air';
params.card = 'Intel HD Graphics 6000'; % find in systems
params.position = 'Hopefully with Miaomiao';

% parameters which can make programming the display easier
params.flipLR = 0;
params.flipUD = 0;

return
