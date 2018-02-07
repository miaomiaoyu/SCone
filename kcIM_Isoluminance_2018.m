% kcIM_Isoluminance_2018

% Display drifting gratings and allows user to adjust theta until stimulus
% appears isoluminant.
% Formula for S cone is [cos(theta)/sqrt(2), cos(theta)/sqrt(2),
% sin(theta)].

% 19-Dec-2017, MY.

clear;
close all;

cd('/Users/miaomiaoyu/GoogleDrive/Matlab_Toolboxes/Projects/Koniocellular');

% - - - - - - - - - - - - - - - - - - - - -
%    Screen Tests: skipped when coding 
% - - - - - - - - - - - - - - - - - - - - -

coding=1;
if coding
    Screen('Preference', 'SkipSyncTests', 1);
    %PsychDebugWindowConfiguration();
end

runningonVP=0;
if runningonVP
    params=displayParamsVP;
else
    params=displayParamsMac;
end

PPD=mmy_calculatePPD(params.distance, params.dimensions(2), ...
    params.numPixels(2));

% - - - - - - - - - - - - - - - - - - - - -
%          Initial Parameters
% - - - - - - - - - - - - - - - - - - - - -

W = what; exptPath = strcat(W.path,'/'); clear W;

coneType = {'L-M', 'S'};
nReps = 3;

prompt = {'Subject ID number: '}; dlgTitle = 'Subject details';
subjID = inputdlg(prompt, dlgTitle, 1);
useSCone = listdlg('PromptString', 'Cone Type:',...
    'SelectionMode', 'single',...
    'ListString', coneType,...
    'ListSize', [160, 40]);

useSCone = useSCone - 1; % if you choose S Cone, useSCone will be 1.
% R-G, useSCone will be 0.

fName = [exptPath, '/kcIM_Isoluminance/S', cell2mat(subjID),...
    '_thetaVals_', coneType{useSCone+1}, '.mat'];

[~, CompName] = system('hostname');
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ...
        && (length(Screen('Screens')) > 1)
    screenNumber = max(Screen('Screens')) - 1;
else
    screenNumber = max(Screen('Screens'));
end

rect = Screen('Rect', screenNumber); % get the screen resolution.
centreX = rect(3) / 2; centreY = rect(4) / 2;
refreshRateHz = Screen('NominalFrameRate', screenNumber); % will be 120 on ViewPixx

% - - - - - - - - - - - - - - - - - - - -
%          Fixation Parameters
% - - - - - - - - - - - - - - - - - - - -
% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
% if you're using a cross instead:

fCross = fixation_cross(2, 10, centreX, centreY);

ringRadius = rect(4) / 2; % outer edge (radius) of the ring: the edge of the screen
ringWidth = ringRadius - PPD / 3; % 1/3 of a degree thick

% Make the ring. It's in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(ringRadius, 1) > 0.5);
imSize = rect(4); % making this the same size as the page

% Define the ring:
xx = (1 - imSize) / 2:(imSize - 1) / 2;
[xx, yy] = meshgrid(xx, xx);
[~, r] = cart2pol(xx, yy);

% make the alpha mask for the ring.
ringAlpha = ((r > ringWidth + 1) & (r < ringRadius - 1));

% - - - - - - - - - - - - - - - - - - - -
%        Stimulus Parameters: DKL
% - - - - - - - - - - - - - - - - - - - -

% Loading in the required parameters:
load('Viewpixx_Processed_cal_data_2_4_2016.mat');

% Stockman/Sharpe cone fundamentals:
load('StockmanSharpe_2deg_cone_fundamentals_1nm.mat');
fundamentals=fundamentals2deg;
gaborSizeDeg = 10;

if useSCone
    maxTheta = deg2rad(135); % 135
    minTheta = deg2rad(45); % 45  0.7854 to 2.3562
    thetaStep = deg2rad(1);
else
    maxTheta = deg2rad(179); % 179 * you want it at about 177?
    minTheta = deg2rad(60); % 90
    thetaStep = deg2rad(2); % 1.5708 to 3.1241
    
    % Alex's is about 1.5708 for LM, 1.6872 for S. 
end

% - - - - - - - - - - - - - - - - - - - -
%    Stimulus Parameters: Gabor Patch
% - - - - - - - - - - - - - - - - - - - -

% Gabor Properties

gaborDimPix = round(PPD * gaborSizeDeg);
sigma = gaborDimPix / 6;
contrast = 1;
aspectRatio = 1.0;
nCycles = 10;
spatialFreq = nCycles / gaborDimPix;
temporalFreq = 20; % not sure if I need to use this if I'm drifting it.
orientation = 0;
backgroundOffset = [0 0 0 0];
disableNorm = 1;
preContrastMultiplier = 1;

% - - - - - - - - - - - - - - - - - - - -
%              Keyboard
% - - - - - - - - - - - - - - - - - - - -

% Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames');

quitCode = KbName('ESCAPE');
pauseCode = KbName('p');
upTheta = KbName('RightArrow');
downTheta = KbName('LeftArrow');
selectTheta = KbName('y');

% - - - - - - - - - - - - - - - - - - - -
%         Initialising Display
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
    if runningonVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        %%%PsychImaging('PrepareConfiguration');
        %%%PsychImaging('AddTask', 'General', 'UseDataPixx');
        
        %%%Datapixx('Open');
        
        %The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        %%%Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        %%%Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        %%%Datapixx('RegWrRd'); % Synchronize Datapixx regiesters to local register cache
        
        %%%PsychImaging('AddTask','General','FloatingPoint42BitIfPossible');
        %%%PsychImaging('AddTask','General','EnableDataPixxM16OutputWithOverlay');
        
        %if Datapixx('IsViewpixx3D') %If it's the Viewpixx3D
        
        % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
        % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
        % The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
        % The gamma values here were obtained following measurements (through the goggles) on
        % the Jaz Spectrometer taken 2/4/2016.
        % We will simply average the left and right eye's values, because they are so similar.
        
        gammaRed = mean([GammaValues(1,1,1), GammaValues(1,1,2)]);
        gammaGreen = mean([GammaValues(2,1,1), GammaValues(2,1,2)]);
        gammaBlue = mean([GammaValues(3,1,1), GammaValues(3,1,2)]);
        
        Screen('LoadNormalizedGamma');
        Screen('LoadNormalizedGammaTable', w, linspace(0,1,256)'* ones(1,3), 0);
        
        % We'll use the average of the right and left gamma values
        PsychColorCorrection('SetEncodingGamma', w, [1/gammaRed, 1/gammaGreen, 1/gammaBlue]);
        
        %%%Datapixx('DisableVideoLcd3D60Hz'); % => According to Daniel B, disabling seems to give less crosstalk, bizarrely!
        %%%subjectData.DisplayType = 'Viewpixx3D'; % set aside the device type for reference
        %%%Datapixx('RegWr');
        
        if Datapixx('IsPropixx') % if it's the Propixx DLP projector
            
            subjectData.DisplayType = 'PROpixx'; % set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram', 0); % set to normal RGB video processing for driving the LEDs & DLP MMDs
            
            Datapixx('RegWrRd'); % seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
            
        end % if Datapixx
        
    end %if runningonVP
    
    % No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    % raise priority level:
    priorityLevel = MaxPriority(w);
    Priority(priorityLevel);
    
    % Set the alpha-blending:
    % We want a linear superposition of the dots should they overlap:
    % Just like the Gabors in GarboriumDemo.m (see there for further info).
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE); % Not sure this is actually working for the color dots)
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    % not sure how this will conflict with the above command
    % about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % ==> definitely need this!!!
    
    % - - - - - - - - - - - - - - - - - - - -
    %         Generating Textures
    % - - - - - - - - - - - - - - - - - - - -
    
    % 1. Gabor
    gaborTexture = CreateProceduralGabor(w, gaborDimPix, gaborDimPix, [], ...
        backgroundOffset, disableNorm, preContrastMultiplier);
    
    % 2. Fixation Ring
    ringMat(:, :, 1) = fixationRing;
    ringMat(:, :, 2) = ringAlpha;
    fRingTexture = Screen('MakeTexture', w, ringMat, [], [], 2);
    
    % - - - - - - - - - - - - - - - - - - - -
    %          Experiment Parameters
    % - - - - - - - - - - - - - - - - - - - -
    
    % Query the screen refresh rate:
    ifi = Screen('GetFlipInterval', w); % Duration between screen flips (1 / refreshRateHz).
    vbl = Screen('Flip', w); % Time it takes to flip.
    
    % - - - - - - - - - - - - - - - - - - - -
    %            Experiment Loop
    % - - - - - - - - - - - - - - - - - - - -
    quitExperiment = 0;
    
    phaseShiftPerFrame = 30;
    phase = 1;
    
    orientationShiftPerFrame = 0.1;
    
    %           Instructions Page
    % - - - - - - - - - - - - - - - - - - - -
    
    commandwindow;
    
    line1='Keep fixation on the cross at the centre of the screen.';
    line2='\nTo QUIT press <ESCAPE>. To PAUSE press <P>, UNPAUSE press <U>.';
    line3='\nUse <LEFT> and <RIGHT> arrows to adjust the gratings until it looks least drifty.';
    line4='\nPress <Y> to confirm your selection.';
    line5='\n\nPress any key to continue.';
    Screen('TextFont', w, 'Courier New');
    Screen('TextSize', w, 30);
    [nx, ny, bbox] = DrawFormattedText(w, [line1 line2 line3 line4 line5], 'center', 'center', [0 0 0]);
    Screen('FrameRect', w, 0, bbox);
    Screen('DrawText', w, '', nx, ny, [255, 0, 0, 0])
    Screen('Flip', w);
    
    KbCheck();
    
    userResponse = 0;
    
    while ~userResponse
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(quitCode)
                sca;
                return;
                disp('*** User chose to quit before experiment has begun. ***');
            else
                userResponse = 1;
            end
        end
    end
    
    %           Stimuli Presentation
    % - - - - - - - - - - - - - - - - - - - -
    % Frames
    f = 0; % this value increases with each iteration, but on a PER-EYE basis only
    missedFrames = 0;
    
    vbl = Screen('Flip', w);
    
    repIndex = 0;
    
    while repIndex < nReps
        
        acceptTheta = 0;
        
        theta = rand(1) * (maxTheta-minTheta) + minTheta; % Randomize a new theta
        disp('this is starting theta:')
        disp(theta)
        % - - - - - - - - - - - - - - - - - - - -
        %            SUPERFAST LOOP
        % - - - - - - - - - - - - - - - - - - - -
        % try to put as little here as possible for efficient
        % processing
        
        lastKeyPressNow = GetSecs; % Last time we accepted a keypress.
        
        % Set the colour:
        if useSCone
            stimLMS.dir = [cos(theta)/sqrt(2), cos(theta)/sqrt(2), sin(theta)]; % should be [0 0 1]
            stimLMS.scale = 0.3;
            
      
        else
            stimLMS.dir = [cos(theta), sin(theta), 0]; % should be [1 -1 0]
            stimLMS.scale = 0.05;
     
        end    % im gna set it to this altho idk if its right. if u replace S with L and L+M with M.
        
        stimRGB = LMS2RGB_Vpixx(stimLMS, fundamentals, resampledSpectra);
        
        while ~acceptTheta
            
            commandwindow;
            
            t = 0;
            
            phase = phase + phaseShiftPerFrame;
            % Here you can, if you want, do something else to make phase either
            % 0 or 180 so that it flickers...
            
            orientation = orientation + orientationShiftPerFrame;
            
            colorMod = stimRGB.dir .* stimRGB.scale*2;
            
            propMat = [phase, spatialFreq, sigma, contrast, aspectRatio, [0 0 0]]';
            
            % - - - - - - - - - - - - - - - - - - - -
            %           Drawing Stimuli
            % - - - - - - - - - - - - - - - - - - - -
            Screen('BlendFunction', w, GL_ONE, GL_ONE);
            
            Screen('DrawTextures', w, gaborTexture, [], [], orientation, [], [], colorMod, [],...
                kPsychDontDoRotation, propMat);
            
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
            
            Screen('DrawTexture', w, fRingTexture);
            Screen('FillRect', w, [0 0 0], fCross);
            
            Screen('BlendFunction', w, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
            Screen('BlendFunction', w, GL_ONE, GL_ZERO);
            
            Screen('DrawingFinished', w);
            
            [vbl , ~ , ~, missed] = Screen('Flip', w, vbl + (ifi * 0.5), [], [], 1); %update display on next refresh (& provide deadline)
            
            if missed > 0
                missedFrames = missedFrames + 1;
            end
            
            % - - - - - - - - - - - - - - - - - - - -
            %          User Response
            % - - - - - - - - - - - - - - - - - - - -
            
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
            
            %if deltaSecs > 0.5 % deltaSecs is time in seconds between each detected keypress.
            %the thing is u can't just put this, because none of the deltaSecs are  > 0.5, and so you
            % kind of never enter this loop
            
            if (keyIsDown)
                
                timeSinceLastPress= GetSecs - lastKeyPressNow;
                
                if (timeSinceLastPress > .2) % This forces 200ms between keypress acceptance
                    
                    lastKeyPressNow = GetSecs;
                    
                    if keyCode(upTheta)
                        
                        theta = theta + thetaStep;
                        %disp(theta);
                        
                        if theta > maxTheta      % cap it at upper limit
                            theta = maxTheta;
                            Beeper(500, 10, 0.05);
                        end
                        % Set the colour:
                        if useSCone
                            stimLMS.dir = [cos(theta)/sqrt(2), cos(theta)/sqrt(2), sin(theta)];
                            stimLMS.scale = 0.30;
                      
                        else
                            stimLMS.dir = [cos(theta), sin(theta), 0];
                            stimLMS.scale = 0.05;
                            
                        end    % im gna set it to this altho idk if its right. if u replace S with L and L+M with M.
                        
                        stimRGB = LMS2RGB_Vpixx(stimLMS, fundamentals, resampledSpectra);
                        
                    elseif keyCode(downTheta) % DECREASE THETA
                        
                        theta = theta - thetaStep;
                        %disp(theta);
                        
                        if theta < minTheta     % Cap it at lower limit 
                            theta = minTheta;
                            Beeper(500, 10, 0.05);
                        end
                        
                        % Set the colour:
                        
                        if useSCone
                            stimLMS.dir = [cos(theta)/sqrt(2), cos(theta)/sqrt(2), sin(theta)];
                            stimLMS.scale = 0.30;
                         
                        else
                            stimLMS.dir = [cos(theta), sin(theta), 0];
                            stimLMS.scale = 0.05;
                            
                        end    % im gna set it to this altho idk if its right. if u replace S with L and L+M with M.
                        
                        stimRGB = LMS2RGB_Vpixx(stimLMS, fundamentals, resampledSpectra);
                        
                    elseif keyCode(selectTheta) % YES THETA
                        
                        repIndex = repIndex + 1;
                        acceptTheta = 1;
                        
                        thetaVals(repIndex) = theta; % It will be randomized again on the next loop
                        
                    elseif keyCode(quitCode)
                        quitExperiment = 1;
                        sca;
                        disp('*** User chose to quit experiment. The existing data has been saved. ***')
                        return
                    end % End check on last keypress time
                end % End check on whether key was pressed
            end
        end
        %end
    end
    
    %%%%------------------------------------------------%%%%
    %               Close down screen:
    %%%%------------------------------------------------%%%%
    %turn off the prioritisation:
    Priority( 0 ); %restore priority
    %Close down the screen:
    Screen('CloseAll');
    
    Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
    
    %Bring back the mouse cursor:
    ShowCursor();
    
catch
    psychrethrow(psychlasterror)
    
    if runningonVP        % close down the ViewPixx or ProPixx
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        %Datapixx('Close'); %closing it here might cause it to crash?
    end
    
    %Close down the screen:
    Screen('CloseAll');
    
    Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
    
    %Bring back the mouse cursor:
    ShowCursor();
    
end % terminating the try/catch statement.

%% Saving out your data

meanTheta = mean(thetaVals);
save(fName, 'meanTheta', 'thetaVals');


%% Sanity check?
% meanTheta for L-M should be around 135-180.
% meanTheta for S cone should be around 45-90.
