%adapted from RM's MID_Dots_test_colour2
%Modification of basic MID_Dots function to get S-cone (& other colour) working!!

%16/4/16
%2nd version, using LMS2RGB_Viewpixx method following failed/aborted use of DKL2RGB.
%Notes:
% The apparent 'hard edge' of the inner annulus was probably only due to having a large dot sigma
% & hence large dot size that meant the two cosine ramps were overlapping % summing, rather than being
% a smooth progression - nothing to worry about as long as the dot sigma stays small!
clear;close all;
PsychDebugWindowConfiguration();
%--------------------------
UseS_cone=1; %Flag for S cone stimulus or L-M cone stimulus (=1 for S cone)
cone_contrast=1; %max for s-cones is 0.748 on the viewpixx
% Which ever one we choose (S or L-M) will have this much contrast
%------------
%Load in Display Params
params=displayParamsVP; %this is the one for viewpixx

[PPD,viewingAngleDeg]=mmy_calculatePPD(params.distance, params.dimensions(2), params.numPixels(2));
%At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)

%If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = true; %Seems better & less flickery when 'true'?

%Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

%Choose the screen: it is usually the max screen no. available.
%Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
%So if we're on that machine, need to -1 from the screen number:
[~, CompName] = system('hostname'); %find out the computer name
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... %normal strcmp not working here, can't figure out why...
        && (length(Screen('Screens'))>1) %and there is more than 1 display connected...
    WhichScreen = max( Screen( 'Screens' ) )-1;
else
    WhichScreen = max( Screen( 'Screens' ) ); %should be right for any other machine!
end

screenRect = Screen('Rect',WhichScreen); %get the screen resolution.
centreX = screenRect(3)/2;
centreY = screenRect(4)/2;
RefreshRate = Screen('NominalFrameRate', WhichScreen);

%jheapcl; %clear the java heap space.

% Define the dot texture, a square-shaped sheet of dots.
%Make the texture the same size as the height of the screen
imsize = screenRect(4);

%set up the raised cosine annular window.
%specify parameters for the annulus:
inrad = PPD * 1;% inner radius of annulus (in pixels), for fixation spot
outrad = PPD * 12/2; %outer radius of annulus (in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = 5; %make it one dot size wide
%This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2; %double the texture size
x0 = (imsize2+1)/2;
y0 = (imsize2+1)/2;
J = ones(imsize2);
for (ii=1:imsize2)
    for (jj=1:imsize2)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > outrad^2)
            J(ii,jj) = 0;
        elseif (r2 < inrad^2)
            J(ii,jj) = 0;
        elseif (r2 > (outrad - cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (inrad + cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end
%%%%-------------------------------%%%%
%       Set up fixation
%%%%-------------------------------%%%%

%Set up the fixation cross or spot:
%This is drawn directly to the screen using Screen('FillRect')
%if you're using a cross instead:
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

%Make the fixation lock ring:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.5;                % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;        % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/4; % 1/4 of a degree thick
ringWidthOuter = ringRadiusOuter - PPD/3; % 1/3 of a degree thick

%Make the ring. It's in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
%Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.
ring_alphaOuter = ((r>ringWidthOuter+1) & (r<ringRadiusOuter-1));

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

%------------------------------------%
% Prepare DKL stimulus

% Need to first load the processed cal data for the Viewpixx:
% It's easier to load it here once rather than every time the RGB conversion runs
% This contains the resampledSpectra.
load('Viewpixx_Processed_cal_data_2_4_2016.mat')

% Load the Stockman/Sharpe 2 deg cone fundamentals:
%load('StockmanSharpe_2deg_cone_fundamentals_1nm.mat')

% Load the Stockman/Sharpe 10 deg cone fundamentals:
load('StockmanSharpe_10deg_cone_fundamentals_1nm.mat')


%Go through each pixel in the dot and assign appropriate RGB values using DKL2RGB
%Compute RGB values for every pixel in LMS cone activation space.
if UseS_cone
    stimLMS.dir = [1 -2 0]; %L-M    the desired LMS cone activation
    stimLMS.scale=.06;
    
    %stimLMS.dir = [ 0 0 1]; %S   the desired LMS cone activation
    %stimLMS.scale=.6;
    
else
    stimLMS.dir = [1 1 1];
        stimLMS.scale=.99;
end

DesiredScale = cone_contrast; %the scale factor (akin to desired cone contrast, I think); max for s-cones is about 0.748 on the viewpixx

%Also plop an alpha channel in to the 4th value of the 3rd dimension.
%Note that these are identical for the 2 polarity dots and describe the transparency values of each
%1=opaque, 0=fully transparent.
% Alpha stu

try %Start a try/catch statement, in case something goes awry with the PTB functions
    
    %----------------------------
    % Set up the screen
    %----------------------------
    
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
    
    
    [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
    UsingVP=0; %take this off when youre actually using vp
    %Initialise the Vpixx device:
    if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        %The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        
        if Datapixx('IsViewpixx3D') %If it's the Viewpixx3D
            
            % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
            % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
            % The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
            % The gamma values here were obtained following measurements (through the goggles) on
            % the Jaz Spectrometer taken 2/4/2016.
            % We will simply average the left and right eye's values, because they are so similar.
            
            gammaRed =   [3];    %mean([GammaValues(1,1,1), GammaValues(1,1,2)]);
            gammaGreen = [3];  %mean([GammaValues(2,1,1), GammaValues(2,1,2)]);
            gammaBlue =  [3]; %mean([GammaValues(3,1,1), GammaValues(3,1,2)]);
            % Screen('LoadNormalizedGamma' ) 
            % We'll use the average of the right and left gamma values
              PsychColorCorrection('SetEncodingGamma', win, [1/gammaRed, 1/gammaGreen, 1/gammaBlue]);
            
            % The following is the gamma used for the achromatic stim in
            % prior versions:
            %K_gamma_left = 2.9712;
            %K_gamma_right = 3.4363;
            % We'll use the average of the right and left gamma values
            %PsychColorCorrection('SetEncodingGamma', win, 1/mean([K_gamma_left, K_gamma_right]));
            
            Datapixx('DisableVideoLcd3D60Hz'); %=> According to Daniel B, disabling seems to give less crosstalk, bizarrely!
            subjectData.DisplayType = 'Viewpixx3D'; %set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx') %if it's the Propixx DLP projector
            
            subjectData.DisplayType = 'PROpixx'; %set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); %set to normal RGB video processing for driving the LEDs & DLP MMDs
          
            Datapixx('RegWrRd'); %seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    end
    %No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    
    %HideCursor;
    %raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    %Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); %in sec
    
    %Set the alpha-blending:
    %We want a linear superposition of the dots should they overlap:
    %Just like the Gabors in GarboriumDemo.m (see there for further info).
    %Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE); %(Not sure this is actually working for the color dots)
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    %not sure how this will conflict with the above command
    %about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %==> definitely need this!!!
    
    
% Convert to LMS stim
   %Do positive polarity dot first:
       
        stimRGB = LMS2RGB_Vpixx (stimLMS, fundamentals10deg, resampledSpectra);
     

% You could really, instead, use the procedural sine wave generation
% stuff from PTB to make these sine waves properly. To flicker - you have
% several options. One is to make two version sof the wave and alternate
% betweeen them (+pi radian). Or you can specify phase into the procedural
% one - that allows you to make drifting gratings etc as well.
[gratingID, gratingRect] = CreateProceduralSineGrating(win, 200, 200,[.5 .5 .5  1 ]);% [, backgroundColorOffset =(0,0,0,0)] [, radius=inf][, contrastPreMultiplicator=1])
colorMod=(stimRGB.dir*stimRGB.scale); 


    %LMSText=Screen('MakeTexture',win,LMSStack);
    
    
    
    % Generate the Inner ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaInner;
    fixationRingTextureInner = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Generate the Outer ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaOuter;
    fixationRingTextureOuter = Screen('MakeTexture',win,ringMat,[],[],2);
    
    f = 0; %this value increases with each iteration, but on a per eye basis only
    missedFrames = 0;
    vbl = Screen('Flip',win); %sync vbl to start time
        % We'll use the average of the right and left gamma values
   % PsychColorCorrection('SetEncodingGamma', win, [1/gammaRed, 1/gammaGreen, 1/gammaBlue]);

    while ~KbCheck
        
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
        % Screen('DrawTexture',win, fixationRingTextureInner);
         Screen('DrawTexture',win, fixationRingTextureOuter);
        %Screen('DrawTexture',win, LMSText );
        phase=(mod(f,10)>5)*180;
        %Screen('FillRect',win,[0 0 0],fixationCross);
        Screen('DrawTexture',win, gratingID,[],[],180,[],[],colorMod*3,[],[],[phase ,.02, .5 , 0]);
        %Draw the fixation cross:
        
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        
        
        Screen('DrawingFinished', win);
        
        [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
        
        f = f+1 %increment counter for next frame
        %keep record of any missed frames:
        if missed > 0
            missedFrames = missedFrames + 1;
        end
    end % End KbCheck
    
    missedFrames
    
    %%%%------------------------------------------------%%%%
    %               Close down screen:
    %%%%------------------------------------------------%%%%
    %turn off the prioritisation:
    Priority( 0 ); %restore priority
    
    if UsingVP        % close down the ViewPixx or ProPixx
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        %Datapixx('Close'); %closing it here might cause it to crash?
    end
    
    %Close down the screen:
    Screen('CloseAll')
    
    Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
    
    %Bring back the mouse cursor:
    ShowCursor();
    
catch
    
    psychrethrow(psychlasterror)
    if UsingVP        % close down the ViewPixx or ProPixx
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        %Datapixx('Close'); %closing it here might cause it to crash?
    end
    
    %Close down the screen:
    Screen('CloseAll')
    
    Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
    
end %End of try/catch statement




