%mmy_FlickerTestColor
%flickers a gabor grating for about 3 seconds then quits automatically
%written by MY, 18-Sep-2017.
% -------------------------------------------------------------------------

%% Clear the workspace
clear; close all;
sca; 

%% Chromatic Conditions
Lum.name='Lum';
Lum.dir=[1 1 1];
Lum.scale=.50;

SCone.name='SCone';
SCone.dir=[0 0 1];
SCone.scale=.60;

LM.name='LM';
LM.dir=[1 -2 0];
LM.scale=.20;

stimLMS={Lum;SCone;LM};

load('Viewpixx_Processed_cal_data_2_4_2016.mat'); 
load('StockmanSharpe_10deg_cone_fundamentals_1nm.mat');

%% Screen Configurations

% Set screenNum, 'max' allows it to connect to any external monitor
screenNum=max(Screen('Screens')); 

% Psych Debugging
%PsychDebugWindowConfiguration();

% Colors
white=WhiteIndex(screenNum);                                     
gray=white/2;

% Window information
[win, winRect] = PsychImaging('OpenWindow', screenNum, gray);

% Flip the screen to clear
Screen('Flip', win);
ifi=Screen('GetFlipInterval', win); % refresh period in seconds

% Parameter(s) for later
dispTimeSec=3;                 % display time
gaborFlickerRateHz=30;                          % flicker time

%% Gabor Properties

gaborDimPix=winRect(4)/4; 
sigma=gaborDimPix/8;

% Basic Gabor features
orientation=0;
contrast=1;
aspectRatio=1.0;
phase=0;
numCycles=10;
freq=numCycles/gaborDimPix; 
backgroundOffset=[0.5 0.5 0.5 0.5];
disableNorm=1;
preContrastMultiplier=0.5;

% Building a procedural gabor texture
gaborTexture=CreateProceduralGabor(win, gaborDimPix, gaborDimPix, [], ...
    backgroundOffset, disableNorm, preContrastMultiplier);

% Building a property matrix to randomise the gabor phases. 
propMat=[phase, freq, sigma, contrast, aspectRatio, 0, 0, 0]';

%% Present the Gabor Patch

% You'll have 3 variables: framesTotal, framesPerStim, and frameCount. What
% you're doing is simple - you're adding framesPerStim to frameCount in
% each loop, then checking if this equals to framesTotal. Once
% frameCount==framesTotal (i.e. you've displayed all your stimuli), the programme quits.


framesTotal=round(dispTimeSec/ifi); % 3s * 60frames/s = 180
framesPerStim=round((1/gaborFlickerRateHz)/ifi); 

startTime=GetSecs; % measures the starting time of session

frameCount=1; % let frame counter start at 0

% Here we can compute a stimulus sequence before we start. It is useful to
% do this so that we can a) save it out later to check it b) it is faster
% c) it can cope with more complex elements like contrast or something.
stimSeq=1:framesTotal;
stimSeq=mod(stimSeq,framesPerStim);
stimSeq=(stimSeq<(framesPerStim/2));

colorMod={[255 0 0],[0 255 0],[0 0 255]};

for i=1:3
    %stimRGB=LMS2RGB_Vpixx(stimLMS{2},fundamentals,resampledSpectra);
    %colorMod=stimRGB.dir*stimRGB.scale;
   for frameCount=1:framesTotal/3
        
        if (stimSeq(frameCount)==1) % if these divide with no remainder
            Screen('DrawTextures', win, gaborTexture, [], [], orientation, [], [], colorMod{i}, [],...
                kPsychDontDoRotation, propMat); % draw the gabor
        end
        
        Screen('Flip', win);
        
        frameCount=frameCount+1;
    end
end

endTime=GetSecs;
totalTime=endTime-startTime

Screen('CloseAll'); 

sca;