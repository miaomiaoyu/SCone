%% KCIM_Experiment_With_Breaks_Tyrion
% - - - - - - - - - - - - - - - - - - - - -

clear; close all;

% !!! ALERTS BEFORE TESTING !!!

% Make sure coding=0; runningonVP=1;
% Display flickering gabor patch and requires user to fixate centrally.
% Central 2 degrees replaced by grey gauss blob.
% Saves out:
% 1) .mat with all conditions and missed frames
% 2) the script itself after running

% 9/10/17, M.Y.
% last edited 01/03/18, M.Y.

% Directories to start with:

thisComputer = computer;

if strcmp(thisComputer, 'MACI64')
    curDir = ('/Users/miaomiaoyu/Documents/GitHub/Koniocellular');
    isoDir=('Users/miaomiaoyu/Documents/GitHub/NeuralOscillations/Koniocellular/KCIM_Isoluminance/S');
    params=displayParamsMac;
    coding=1; runningonVP=0;
else
    curDir = ('/Users/tyrion/Documents/MATLAB/Miaomiao/KoniocellularEEG');
    isoDir=('/Users/tyrion/Documents/MATLAB/Miaomiao/KoniocellularEEG/KCIM_Isoluminance/S');
    params=displayParamsVP;
    coding=0; runningonVP=1;
end

if coding
    Screen('Preference', 'SkipSyncTests', 1);
    nReps = 1;
    blockDurationSec = 8;
    isi = 0.5;
    nBreakBlock=2;
else
    nReps = 5;
    blockDurationSec = 12;
    isi = 1;
    nBreakBlock=3;
end

PPD=mmy_Calculate_PPD(params.distance, params.dimensions(2), ...
    params.numPixels(2));

% Directories to save data out to

if ~exist ([pwd, '/KCIM_Data'], 'file')
    mkdir ([pwd, '/KCIM_Data']);   % Make sure there's a place you can save your data
end

% Prompt subject ID input by user

prompt=[{'Subject ID Number: '}, {'Subject Initials: '}]; dlgTitle='Subject Details';
options.Resize='off';
options.WindowStyle='normal';
options.Interpreter='none';
subjInfo=inputdlg(prompt, dlgTitle, 1, {'', ''}, options);
fName=[pwd, '/KCIM_Data/kcIM_S', subjInfo{1}, '_', subjInfo{2} ...    % make this into Practice_Data for when you're coding
    '_', datestr(now,'dd-mm-yyyy_HHMMSS')];

[~, CompName] = system('hostname');
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ...
        && (length(Screen('Screens')) > 1)
    screenNumber = max(Screen('Screens')) - 1;
else
    screenNumber = max(Screen('Screens'));
end

rect=Screen('Rect', screenNumber); % get the screen resolution.
centreX=rect(3)/2;
centreY=rect(4)/2;
refreshRateHz=Screen('NominalFrameRate', screenNumber); % 120 on ViewPixx

%%    Stimulus Parameters: Gabor Patch
% - - - - - - - - - - - - - - - - - - - -

G.contrast=1;
G.aspectRatio=1;
G.orientation=90;
G.backgroundOffset=[0 0 0 0];
G.disableNorm=1;
G.contrastPreMultiplicator=1;
G.gaborDeg=35;
G.gaborDimPix=round(PPD * G.gaborDeg);
G.nCycles=G.gaborDimPix / 30;
G.gaborSigma=G.gaborDimPix/7;
G.gaborSf=G.nCycles/G.gaborDimPix;

% Load the Stockman/Sharpe 10 deg cone fundamentals:
load('StockmanSharpe_10deg_cone_fundamentals_1nm.mat');

gray=GrayIndex(screenNumber);
sizeOfSquare=2 * PPD; % pretty sure this corresponds to size in pixels...
sizeOfBlob=100; % this is arbitrary, i think. it does make the blob smaller, but not by much.
transLayer=2;
blobSigma=700;
[x,y]=meshgrid(-sizeOfSquare:sizeOfSquare, -sizeOfSquare:sizeOfSquare);
greyBlob=uint8(ones(2*sizeOfSquare+1, 2*sizeOfSquare+1, transLayer) * gray);
size(greyBlob);

% Layer 2 (Transparency aka Alpha) is filled with gaussian transparency
% mask.
xsd=sizeOfSquare/2.0;
ysd=sizeOfSquare/2.0;
greyBlob(:,:,transLayer)=uint8(round(-sizeOfBlob + exp(-((x/xsd).^2)-((y/ysd).^2))* blobSigma));


%%          Fixation Parameters
% - - - - - - - - - - - - - - - - - - - -
% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
% if you're using a cross instead:
fCross = fixation_cross(2, 10, centreX, centreY);

crossColorMat={[1 0 0], [0 0 1]};
crossColor = round(rand(1));

%%        Stimulus Parameters: DKL
% - - - - - - - - - - - - - - - - - - - -

% Need to first load the processed cal data for the Viewpixx:
% It's easier to load it here once rather than every time the RGB conversion runs
% This contains resampledSpectra which is needed in defining colorMod.
load('Viewpixx_Processed_cal_data_2_4_2016.mat');

thetaSConeFile=fullfile([isoDir, subjInfo{1}, '_thetaVals_S.mat']);
thetaLMFile=fullfile([isoDir, subjInfo{1}, '_thetaVals_L-M.mat']);

if ~exist( thetaSConeFile, 'file' ) || ~exist ( thetaLMFile, 'file' ) % Check that these files exist.
    disp('ERROR: Cannot locate subject''s isoluminance theta value.');
    sca
    return
    
else
    load(thetaSConeFile);
    meanThetaSCone = meanTheta;
    
    load(thetaLMFile);
    meanThetaLM = meanTheta;
end

% --- Chromatic Levels ---

% Note: contrast levels of Y3 proj was 30, 30, 5.

lum.name = 'Lum';
lum.dir = [1 1 1];
lum.scale = .99;
lum.trig = 5; % 7, 3, 11

sCone.name = 'SCone';
sCone.dir = [cos(meanThetaSCone)/sqrt(2), cos(meanThetaSCone)/sqrt(2), sin(meanThetaSCone)]; % [0 0 1]
sCone.scale = .729;
sCone.trig = 7;

lm.name = 'LM';
lm.dir = [cos(meanThetaLM), sin(meanThetaLM), 0];  % [1 -2 0]
lm.scale = .232;
lm.trig = 11;


stimLMS = {lum; sCone; lm};
nColor = length(stimLMS);

% --- Temporal Levels ---

stimFreq = [5, 12, 16]; % flicker rates of gabor patch

nFreq = length(stimFreq);

%%          Create Conditions
% - - - - - - - - - - - - - - - - - - - -

condIndex=1;

for i=1:nColor
    for ii=1:nFreq
        stimConditions(condIndex,:)={stimLMS{i}, stimFreq(ii)};
        condIndex=condIndex+1;
    end
end

nCond = length(stimConditions);

%%              Keyboard
% - - - - - - - - - - - - - - - - - - - -

% Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames');
quitCode = KbName('ESCAPE');
pauseCode = KbName('p');
attentionCode = KbName('space');

%%              Triggers
% - - - - - - - - - - - - - - - - - - - -

breakStart=22;
breakEnd=29;
trigSec=1;
crossChangedColour = 0;

%%         Initialising Display
% - - - - - - - - - - - - - - - - - - - -

try % Start a try/catch statement, in case something goes awry with the PTB functions
    % Set up the screen
    % initialization of the display
    AssertOpenGL;
    
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % required for gamma correction through the PsychImaging pipeline:
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    %Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    
    % Open an on screen (grey) window and configure the imaging pipeline
    % Info about the 'blueline' mechanism for synching to the 3D glasses:
    % There seems to be a blueline generation bug on some OpenGL systems.
    % SetStereoBlueLineSyncParameters(windowPtr, windowRect(4)) corrects the
    % bug on some systems, but breaks on other systems.
    % We'll just disable automatic blueline, and manually draw our own bluelines!
    
    [w, wRect] = PsychImaging('OpenWindow', screenNumber, 0.5);
    [screenXPix, screenYPix] = Screen('WindowSize', w);
    
    %Initialise the Vpixx device:
    
    if runningonVP % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        
        Datapixx('Open');
        
        % The following commands are included in demos that apparently work for both the
        % Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight'); % optionally, turn it off first, in case the refresh rate has
        % changed since startup
        Datapixx('EnableVideoScanningBacklight'); % Only required if a VIEWPixx.
        Datapixx('RegWrRd'); % Synchronize Datapixx regiesters to local register cache
        
        PsychImaging('AddTask','General','FloatingPoint42BitIfPossible');
        PsychImaging('AddTask','General','EnableDataPixxM16OutputWithOverlay');
        
        Datapixx('DisableVideoLcd3D60Hz'); % => According to Daniel B, disabling seems to give less crosstalk, bizarrely!
        subjectData.DisplayType = 'Viewpixx3D'; % set aside the device type for reference
        Datapixx('RegWr');
        
        % No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    end
    
    if Datapixx('IsPropixx') % if it's the Propixx DLP projector
        
        subjectData.DisplayType = 'PROpixx'; % set aside the device type for reference
        Datapixx('SetPropixxDlpSequenceProgram', 0); % set to normal RGB video processing for driving the
        % LEDs & DLP MMDs
        Datapixx('RegWrRd'); % seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        
    end % if Datapixx
    
    % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
    % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
    % The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
    % The gamma values here were obtained following measurements (through the goggles) on
    % the Jaz Spectrometer taken 2/4/2016.
    % We will simply average the left and right eye's values, because they are so similar.
    
    gammaRed=mean([GammaValues(1,1,1), GammaValues(1,1,2)]);
    gammaGreen=mean([GammaValues(2,1,1), GammaValues(2,1,2)]);
    gammaBlue=mean([GammaValues(3,1,1), GammaValues(3,1,2)]);
    
    Screen('LoadNormalizedGamma');
    Screen('LoadNormalizedGammaTable', w, linspace(0, 1, 256)'* ones(1, 3), 0);
    
    % We'll use the average of the right and left gamma values
    PsychColorCorrection('SetEncodingGamma', w, [1/gammaRed, 1/gammaGreen, 1/gammaBlue]);
    
    % Raise priority level:
    priorityLevel = MaxPriority(w);
    Priority(priorityLevel);
    
    %%            Alpha Blending
    % - - - - - - - - - - - - - - - - - - - -
    
    % Set the alpha-blending:
    % We want a linear superposition of the dots should they overlap:
    % Just like the Gabors in GarboriumDemo.m (see there for further info).
    
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE); % Not sure this is actually working for the color dots)
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    % not sure how this will conflict with the above command
    % about the linear superposition of dots... but it doesn't seem to cause problems
    
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % ==> definitely need this!!!
    
    %%         Generating Textures
    % - - - - - - - - - - - - - - - - - - - -
    
    gaborTexture=CreateProceduralGabor(w, G.gaborDimPix, G.gaborDimPix, [],...
        G.backgroundOffset, G.disableNorm, G.contrastPreMultiplicator);
    
    % Build a single transparency mask texture
    greyBlobTexture=Screen('MakeTexture', w, greyBlob);
    
    %%          Experiment Parameters
    % - - - - - - - - - - - - - - - - - - - -
    
    quitExperiment=0;
    pauseExperiment=0;
    globalIndex=1;
    
    % Query the screen refresh rate:
    ifi = Screen('GetFlipInterval', w); % Duration between screen flips (1 / refreshRateHz).
    vbl = Screen('Flip', w); % Time it takes to flip.
    
    %%            Experiment Loop
    % - - - - - - - - - - - - - - - - - - - -
    
    lastKeyPressNow = GetSecs;
    
    while quitExperiment<1 && pauseExperiment<1
        
        commandwindow; % this might only work in Matlab 2016 and later.
        
        % Create an 'Instructions' page.
        
        instrText = 'Keep fixation on the cross. \n\n There will be a short break every 10 minutes. \n\nPress <ESCAPE> to QUIT. \n\n Press <ANY KEY> to start.';

        [nx, ny, bbox] = DrawFormattedText(w, instrText, 'center', 'center', [0 0 0]);
        
        Screen('FrameRect', w, 0, bbox);
        
        Screen('DrawText', w, '', nx, ny, [255, 0, 0, 0])
        
        Screen('Flip',w);
        
        KbWait;
        
        % -------- F I X A T I O N  C R O S S  --------
        %
        %         for secondChunk = 1:5
        %             crossColor = round(rand(1));
        %             Screen('FillRect', w, crossColorMat{crossColor+1}, fCross);
        %         end
        %
        
        % -------- F I X A T I O N  C R O S S  --------
        
        for breakIndex=1:nBreakBlock % for each 5 reps you get a break
            
            for repIndex=1:nReps
                
                displayIndex=randperm(nCond); % Gives a random order (since we can't randperm cells directly)
                
                for condIndex=1:nCond
                    
                    missedFrames=0;
                    % Two things get set here : The color (in LMS) and the temporal
                    % frequency (from the array gaborFlickerRateHz).
                    % For each presentation of the stim, you want to work out how
                    % many frames there are. Then compute a sine wave that has that
                    % many frames in it and where the frequency is changing at x HZ
                    % (set by the gaborFlickerRateHz).
                    
                    % First do the color
                    stimRGB = LMS2RGB_Vpixx(stimConditions{displayIndex(condIndex), 1}, fundamentals10deg, resampledSpectra);
                    colorMod = stimRGB.dir .* stimRGB.scale;
                    
                    % Now do the stim seq
                    framesPerStim=round(blockDurationSec / ifi); % 6s * 120frames/s = 720 frames
                    thisFreq=stimConditions{displayIndex(condIndex), 2};
                    cyclesPerBlock=thisFreq * blockDurationSec;
                    modAngle=linspace(0, 2 * pi * cyclesPerBlock, framesPerStim);
                    
                    modWave=sin(modAngle); % Contrast Reversing Sine wave
                    modWave=sign(modWave); % CR Square wave (looks clearer) - ViewPixx prefers this?
                    modWave(modWave<0)=0; % On/Off
                    
                    % Create a unique Condition Code
                    condCode=mmy_Create_Cond_Code(stimConditions { displayIndex(condIndex) }.trig, ...
                        stimConditions { displayIndex(condIndex), 2 } );

                    gaborPhase = rand(1) * 360; % Randomize the phase to avoid cone-level adaptation
                    % if the same part of the screen is always
                    % green/blue/black whatever you'll slowly
                    % come to expect it/adapt to it.
                    
                    frameNo=0;
                    
                    stimStartTime=GetSecs;
                    
                    if runningonVP   % TRIGGER!!! CONDCODE tiny bit BEFORE STIMULUS GOES UP
                        if frameIndex == 0
                            Datapixx('SetDoutValues', transformindex(condCode));
                            Datapixx('RegWrRd');
                        end
                    end
                    
                    % ****** VERY FAST LOOP GOES HERE *********
                    
%                     nCrossChange =4;
%                     
%                     for crossChange = 1:nCrossChange
%                         for frame = 1:405 
%                             
%                             %Screen('BlendFunction', w, GL_ONE, GL_ONE);
%                             crossColor = round(rand(1));
%                             Screen('FillRect', w, crossColorMat{crossColor+1}, fCross);
%                             Screen('Flip', w);
%                         end
%                     end
                    i = 1;
                    
                    for frameIndex = 1:framesPerStim % for each block there are this many frames
                        
                        %      KbCheck for quit and pause
                        % - - - - - - - - - - - - - - - - - - - -
                        
                        [keyIsDown, secs, keyCode] = KbCheck; % keep checking for key presses
                        
                        if keyIsDown
                            
                            timeSinceLastPress = GetSecs - lastKeyPressNow;
                            
                            if (timeSinceLastPress > .5) % this forces 500ms between key press acceptance.
                                
                                lastKeyPressNow = GetSecs;
                            
                            if keyCode(pauseCode)
                                
                                pause('on');
                                
                                if runningonVP
                                    Datapixx('SetDoutValues', transformindex(trigPause));
                                    Datapixx('RegWrRd');
                                end
                                
                                pause;
                                
                                % you can't press 'p' again, it'll just keep
                                % pausing. you have to feed it a different value.
                                % Get user to press 'u' instead.
                                if runningonVP
                                    Datapixx('SetDoutValues', transformindex(trigUnpause));
                                    Datapixx('RegWrRd');
                                end
                                
                            elseif keyCode(quitCode)
                                
                                quitExperiment=1;
                                
                                stimEndTime=GetSecs;
                                stimTime=stimEndTime-stimStartTime;
                                
                                missedFrames=0;
                                sca;
                                fprintf('*** User chose to quit experiment. The existing data has been saved. ***');
                                
                                %%        Saving the data out
                                % - - - - - - - - - - - - - - - - - - - -
                                
                                KCIM_Data(1,:) = {'Global Index', 'Rep Index', 'Cond Index', 'Display Index',...
                                    'Color', 'Freq', 'Phase', 'ISI Start', 'ISI End', 'ISI Duration',...
                                    'Stim Start', 'Stim End', 'Stim Time',...
                                    'Missed Frames', 'Cond Code'};
                                
                                KCIM_Data(globalIndex+1, :) = {globalIndex, repIndex, condIndex, displayIndex(condIndex),...
                                    stimConditions{displayIndex(condIndex), 1}.name,...
                                    stimConditions{displayIndex(condIndex), 2},...
                                    gaborPhase, isiStartTime, isiEndTime, isiTime,...
                                    stimStartTime, stimEndTime, stimTime,...
                                    missedFrames, condCode};
                                
                                save(fName,'KCIM_Data');
                                
                                globalIndex=globalIndex+1;
                                
                                return
                                
                            elseif keyCode(attentionCode)
                                crossChangedColour = GetSecs;
                                
                            end
                            
                            attentionKeyPress(i) = crossChangedColour;
                            
                            i = i+1;
                            
                            end
                        
                        end % End check on keycode being down
                        
                        if runningonVP
                            
                            if (mod(frameIndex, refreshRateHz) == 0) % confirmed: sends out trigger '1' at every second interval.
                                Datapixx('SetDoutValues', transformindex(trigSec));
                                Datapixx('RegWrRd');
                            end
                            
                        end
                        
                        %%           Drawing Stimuli
                        % - - - - - - - - - - - - - - - - - - - -
                        
                        if frameIndex == 90
                            pickAColor = true;
                        else
                            pickAColor = false;
                        end
                        
                        if pickAColor
                            crossColor = round(rand(1));
                        end
                                
                        
                        gaborPropMat = [gaborPhase, G.gaborSf, G.gaborSigma, modWave(frameIndex), ...
                            G.aspectRatio, [0 0 0]]';
                        
                        % Drawing Gabor Patch
                        
                        Screen('BlendFunction', w, GL_ONE, GL_ONE);
                        
                        Screen('DrawTexture', w, gaborTexture, [], [], G.orientation, [], [], colorMod, [],...
                            kPsychDontDoRotation, gaborPropMat);
                        
                        % Drawing Grey Blob
                        
                        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        
                        Screen('DrawTexture', w, greyBlobTexture);
                        
                        % Drawing Fixation Cross
                        
                        
                        Screen('FillRect', w, crossColorMat{crossColor+1}, fCross);
                        %Screen('Flip', w);
                        
                        Screen('DrawingFinished', w);
                        
                        % Update display on next refresh (& provide deadline)
                        [vbl , ~ , ~, missed] = Screen('Flip', w, 0, [], [], 1);
                        
                        
                        if missed>0
                            missedFrames=missedFrames+1;
                            missedTime(missedFrames)=GetSecs; % You should also log the start time. And I'm not sure that now is a good clock. PTB prob has a better one.
                        end
                        
                        if runningonVP
                            
                            if (mod(frameIndex, refreshRateHz) == 0) % confirmed: sends out trigger '0' after screen flip ponce a second
                                Datapixx('SetDoutValues', 0);
                                Datapixx('RegWrRd');
                            end
                            
                        end
                        
                    end % End of 1 stimulus
                    
                    % ******** END OF FAST LOOP *********
                    
                    stimEndTime=GetSecs;
                    stimTime=stimEndTime-stimStartTime;
                    
                    %%                 ISI
                    % - - - - - - - - - - - - - - - - - - - -
                    isiStartTime=GetSecs;
                    
                    framesPerISI=round(isi/ifi);
                    
                    for f=1:framesPerISI % nFrames in an ISI
                        
                        Screen('FillRect', w, [0.5 0.5 0.5]);
                        Screen('Flip', w);
                        
                        if runningonVP
                            if f < 3
                                Datapixx('SetDoutValues', transformindex(4)); % isi
                                Datapixx('RegWrRd');
                            end
                        end
                        
                    end
                    
                    isiEndTime = GetSecs;
                    isiTime = isiEndTime - isiStartTime;
                    
                    %%        Saving the data out
                    % - - - - - - - - - - - - - - - - - - - -
                    
                    KCIM_Data(1,:) = {'Global Index', 'Rep Index', 'Cond Index', 'Display Index',...
                        'Color', 'Freq', 'Phase', 'ISI Start', 'ISI End', 'ISI Duration',...
                        'Stim Start', 'Stim End', 'Stim Time',...
                        'Missed Frames', 'Cond Code'};
                    
                    KCIM_Data(globalIndex+1,:) = {globalIndex, repIndex, condIndex, displayIndex(condIndex),...
                        stimConditions{displayIndex(condIndex), 1}.name,...
                        stimConditions{displayIndex(condIndex), 2},...
                        gaborPhase, isiStartTime, isiEndTime, isiTime,...
                        stimStartTime, stimEndTime, stimTime,...
                        missedFrames, condCode};
                    
                    globalIndex=globalIndex+1;
                    
                end % Every 1 rep of 9 conds 
                
                fprintf('\nTotal number of missed frames: %g', missedFrames);
                
            end % Every 5 reps of 9 conds
            
             if breakIndex < nBreakBlock % As long as it's a break before the last set
                
                breakText = 'You may now take a break \n\n The next block will take approximately 10 minutes. \n\n Press <ANY KEY> to continue.';
                
                [nx, ny, bbox] = DrawFormattedText(w, breakText, 'center', 'center', [0 0 0]);
                
                Screen('FrameRect', w, 0, bbox);
                
                Screen('DrawText', w, '', nx, ny, [255, 0, 0, 0])
                
                if runningonVP
                
                Datapixx('SetDoutValues', transformindex(breakStart)); % Trig to indicate when the break starts
                Datapixx('RegWrRd');
                
                end
                
                Screen('Flip',w);
                
                KbWait;
                
                if runningonVP
                    
                Datapixx('SetDoutValues', transformindex(breakEnd)); % Trig to indicate when the break finishes
                Datapixx('RegWrRd');
                
                end
            else
                continue
            end
            
        end
        
        Screen('Flip', w);
        
        save(fName,'KCIM_Data'); % Save out the data
        
        %% Turn off everything
        
        %turn off the prioritisation:
        Priority( 0 ); %restore priority
        
        %Close down the screen:
        Screen('CloseAll');
        
        Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
        
        %Bring back the mouse cursor:
        ShowCursor();
        
    end % While loop
    
catch
    
    psychrethrow(psychlasterror)
    
    if runningonVP
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        % Datapixx('Close'); %closing it here might cause it to crash?
    end
    
    Screen('CloseAll');
    
    Datapixx('Close');
    
    
end

sca;

