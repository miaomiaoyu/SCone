
% Clear the workspace and the screen
sca;
close all;
clearvars;

% Setup PTB with some default values
PsychDefaultSetup(2);

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;

% Skip sync tests for demo purposes only
Screen('Preference', 'SkipSyncTests', 2);

% Open the screen
[w, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
    [], [],  kPsychNeed32BPCFloat);

%--------------------
% Gabor information
%--------------------

% Dimension of the region where will draw the Gabor in pixels
gaborDimPix = windowRect(4) / 2;

% Sigma of Gaussian
sigma = 30;%gaborDimPix / 7;

% Obvious Parameters
orientation = 0;
contrast = 1;
aspectRatio = 1;
phase = 1;

% Spatial Frequency (Cycles Per Pixel)
% One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
numCycles = 5;
freq = numCycles / gaborDimPix;

% Build a procedural gabor texture (Note: to get a "standard" Gabor patch
% we set a grey background offset, disable normalisation, and set a
% pre-contrast multiplier of 0.5.
% For full details see:
% https://groups.yahoo.com/neo/groups/psychtoolbox/conversations/topics/9174
backgroundOffset = [0.5 0.5 0.5 0.5];
disableNorm = 0;
preContrastMultiplier = 0.5;
%gabortex = CreateProceduralGabor(window, gaborDimPix, gaborDimPix, [],...
% backgroundOffset, disableNorm, preContrastMultiplier);

gabortex = CreateProceduralGaussBlob(w, gaborDimPix, gaborDimPix, backgroundOffset,...
    disableNorm, preContrastMultiplier);

% Randomise the phase of the Gabors and make a properties matrix.
propertiesMat = [contrast,sigma, aspectRatio, 0];


%------------------------------------------
%    Draw stuff - button press to exit
%------------------------------------------

% If you've drawn a blue blob (first stimulus), followed by a red blob
% (second stimulus).

% If you put a 'BlendFunction' after another, it overrides it. 

% % Screen('BlendFunction', w, GL_ONE, GL_ONE); % puts up white square
% Screen('BlendFunction', w, GL_ZERO, GL_ONE) % puts up grey square
% Screen('BlendFunction', w, GL_ZERO, GL_ZERO) % puts up black square
% Screen('BlendFunction', w, GL_ONE, GL_ZERO) % red blob (second stimulus)
%Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %
% purplish blob but more towards red...
%Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE) % white square

%Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA); % very
%nicely blended purple blob

Screen('BlendFunction', w, GL_ONE, GL_ONE_MINUS_DST_ALPHA); % very
%nicely blended purple blob

% Draw the Gabor. By default PTB will draw this in the center of the screen
% for us.

Screen('DrawTextures', w, gabortex, [], [], orientation, [], [], [0 0 255], [],...
    kPsychDontDoRotation, propertiesMat');

Screen('DrawTextures', w, gabortex, [], [], orientation, [], [], [255 0 0], [],...
    kPsychDontDoRotation, propertiesMat');

% Flip to the screen
Screen('Flip', w);

% Wait for a button press to exit
KbWait;

% Clear screen
sca;

% Published with MATLAB® R2015b