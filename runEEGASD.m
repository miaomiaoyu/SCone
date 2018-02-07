function runEEGASD

% display gratings at different contrasts for SSVEP
% also include mask conditions
% central detection/discrimination task
% each participant to run FOUR blocks
% DHB 31/12/14

clear
close all;
W = what; E.exptpath = strcat(W.path,'/'); clear W;

[subj, blocksdone, sno] = runsubjectroutine(E.exptpath);
E.subj = subj;
E.blockstorun = 4-blocksdone;

useVP = 1;

E.nconds = 1;
E.nlevels = 14;
E.ntrials = 8;
E.ntrialsperblock = 28;

E.nrepetitions = E.ntrials/(E.ntrialsperblock/(sum(E.nlevels)));


        E.targetlevelsdB = [0:6:36 0:6:36];
        E.targetlevelsC = d_dBtoPercent(E.targetlevelsdB)./100;
        E.targetlevelsC([1 8]) = 0;
        E.masklevelsdB(1:7) = 0;
        E.masklevelsdB(8:14) = 30;
        E.masklevelsC = d_dBtoPercent(E.masklevelsdB)./100;
        E.masklevelsC(1:7) = 0;
        E.maskori = 90;
        E.maskSF = 0.5;

E.targetphase = 90;
E.maskphase = 0;

DP.ntrialspercond = E.nrepetitions*300;
DP.noiselevelsdB = 6;
DP.noiselevelsC = d_dBtoPercent(DP.noiselevelsdB)./100;
DP.ncontrasts = length(DP.noiselevelsdB);
DP.targetlevelsdB = [-18 6];
DP.targetlevelsC = d_dBtoPercent(DP.targetlevelsdB)/100;
DP.targetlevelsC(1) = 0;
DP.nlevels = length(DP.targetlevelsC);
DP.ntrialsperblock = 600;     % must be divisible by E.nlevels
DP.nrepetitions = DP.ntrialspercond/(DP.ntrialsperblock/DP.nlevels);
DP.ISI = 0.4;
DP.ITI = 0.5;
DP.duration = 0.5;
DP.lookaheadtrials = 5;

KbName('UnifyKeyNames');
clear PsychHID;
if ~useVP
    Screen('Preference', 'SkipSyncTests', 1);
end

ST.npixelsperdegree = 36;       % at 57cm
ST.SF = 0.5;
ST.gratingsize = 110;

ST.ncycles = ST.SF*ST.gratingsize/ST.npixelsperdegree;
ST.maskcycles = E.maskSF.*ST.gratingsize./ST.npixelsperdegree;

ST.duration = 11;
ST.ITI = 3;         % always leave 3 seconds between trials
ST.TFtarget = 7;
ST.TFmask = 5;
ST.TFbehavioural = 6;

window = make_soft_window(ST.gratingsize,ST.gratingsize,0.8);
dpstim = mkgrating(ST.gratingsize, ST.ncycles, 90, E.targetphase, 1) .* window;

fname = strcat(E.exptpath, 'Results/', E.subj, 'trialorderMSc15.mat');
if exist(fname)
    load(fname);
    save(strcat(E.exptpath, 'Results/', E.subj, 'trialorderMSc15old.mat'),'R','E','DP');
else
    R.subj = E.subj;
    R.resps = zeros(E.nconds,max(E.nlevels),E.ntrials);
    R.tno = zeros(E.nconds,max(E.nlevels));
    
    R.orientations = 360.*(rand(E.nconds,max(E.nlevels),E.ntrials));
    
    R.triallistordered(1,:) = ceil((1:(E.nconds.*E.nlevels))./E.nlevels);
    R.triallistordered(2,:) = mod((1:(E.nconds.*E.nlevels))-1,E.nlevels)+1;
    
    R.blocktrialsordered = repmat(1:(E.nconds*E.nlevels),[1 E.ntrialsperblock/(E.nconds*E.nlevels)]);
    
    for n = 1:E.nrepetitions
        R.trialorder(n,:) = randperm(E.ntrialsperblock);
    end
    R.currentrep = 0;
    
    % now set up the double pass conditions for the behavioural task
    DP.ntrials = zeros(size(DP.targetlevelsC));
    DP.resps = zeros(DP.nlevels,DP.ntrialspercond,2);      % last index is the pass
    DP.tno = zeros(DP.nlevels,2);
    
    DP.seeds = randn(DP.nlevels,DP.ntrialspercond,2);        % last index is the trial interval
    DP.testintervals = ceil(rand(DP.nlevels,DP.ntrialspercond,2)*2);
    
    for n = 1:DP.nrepetitions
        DP.blockorder((n-1)*DP.ncontrasts+1:n*DP.ncontrasts) = randperm(DP.ncontrasts);
    end
    DP.currentblockindex(1:2) = 0;
    
    for n = 1:(DP.nlevels*DP.ntrialspercond)
        currenttrial = 0;
        while currenttrial==0				% select current condition
            trialcond = ceil(rand*DP.nlevels);	% choose condition for trial
            
            if DP.ntrials(trialcond)<DP.ntrialspercond   % keep equal n of trials per level
                currenttrial = trialcond;
            end
        end
        DP.trialorder(n) = currenttrial;
        DP.ntrials(trialcond) = DP.ntrials(trialcond) + 1;
    end
    
    save(fname, 'R', 'E', 'DP');
end

WaitSecs(0.01);                % important to load in the MEX file before the expt starts
GetSecs;
InitializePsychSound;
tr = PsychPortAudio('Open',[],[],[],[],1);
PsychPortAudio('FillBuffer', tr, ones(1,220));
tc = PsychPortAudio('Open',[],[],[],[],1);
PsychPortAudio('FillBuffer', tc, MakeBeep(440,0.05,44100).*0.5);

try                 % start the 'try/catch' loop
    
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    PsychGPUControl('SetDitheringEnabled', 0);
    
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    if useVP        % using a ViewPixx or ProPixx
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        Datapixx('EnableVideoScanningBacklight');       % Only required if a VIEWPixx.
        Datapixx('RegWr');
        
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
        [w, winRect] = PsychImaging('OpenWindow', screenNumber, 0, [], [], [], 0, 0, kPsychNeedFastBackingStore);
        Screen('LoadNormalizedGammaTable', w, linspace(0,1,256)'*ones(1,3), 0); % THIS IS THE IMPORTANT THING TO DO, NOTE THE LAST ARGUMENT IS 0.
        
        HideCursor;
        ST.greylevel = doimagegamma(0.5);
    else
        rect = [1 1 1024 1024];
        [w, winRect] = Screen('OpenWindow',screenNumber,128,rect);
        ST.greylevel = 128;
    end
    
    
    [width, height] = Screen('WindowSize', w);
    ifi=Screen('GetFlipInterval', w);
    ifims = ifi*1000;
    
    Screen('TextSize', w, 18);
    
    Screen('FillRect', w, ST.greylevel);
    drawfixation(w, width, height, 0);
    Screen('Flip', w);
    
    ST.nframes = round(1000/ifims);       % no of frames in 1 sec
    
    p = 270*pi/180;         % force waveform to start at a minima
    targetframes = ST.nframes;
    targetwaveform = sin(p + 2 .* ST.TFtarget .* (ifims:ifims:targetframes*ifims) .* pi./1000);
    targetwaveform = (targetwaveform + 1)./2;
    maskframes = ST.nframes;
    maskwaveform = sin(p + 2 .* ST.TFmask .* (ifims:ifims:maskframes*ifims) .* pi./1000);
    maskwaveform = (maskwaveform + 1)./2;
    
    dpframes = 3*round(1000/ifims)/ST.TFbehavioural;
    dpwaveform = sin(p + 2 .* ST.TFbehavioural .* (ifims:ifims:dpframes*ifims) .* pi./1000);
    dpwaveform = (dpwaveform + 1)./2;
    
    DP.nframesinISI = round(1000*DP.ISI/ifims);
    DP.nframesinISI = ceil(DP.nframesinISI/2)*2;        % make sure this is an even number of frames
    DP.nframesinITI = round(1000*DP.ITI/ifims);
    DP.nframesinITI = ceil(DP.nframesinITI/2)*2;        % make sure this is an even number of frames
    
    comp = zeros(ST.gratingsize);
    comp = (1+comp)/2;
    
    if useVP
        comp = doimagegamma(comp);
    end
    nulltexture = Screen('MakeTexture', w, comp, [], [], 2);
    
    r1 = [1 1 ST.gratingsize ST.gratingsize];
    rectcount = 0;
    for x = 1:17
        for y = 1:9
            rectcount = rectcount + 1;
            offsetx = ST.gratingsize*(x-9);
            offsety = ST.gratingsize*(y-5);
            if x==9 && y==5
                rectcount = rectcount-1;
            else
                destRects(:,rectcount) = CenterRectOnPoint(r1, width*0.5+offsetx, height*0.5+offsety);
            end
        end
    end
    dpRect = CenterRectOnPoint(r1, width*0.5, height*0.5);
    
    alltriggertimes = [];
    allframetimes = [];
    
    breakcode = 0;
    currentblock = 0;
    
    while currentblock < E.blockstorun
        
        currentblock = currentblock + 1;
        R.currentrep = R.currentrep + 1;
        
        R.reportedifi(R.currentrep) = ifi;
        
        if currentblock>1
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, width, height,0)
            textstring = strcat('Block', num2str(currentblock-1),' complete');
            Screen('DrawText', w, textstring, (width/2)-200, (height/2)-400, [0 0 0]);
            Screen('Flip', w);
            WaitSecs(10);        % force subjects to wait between blocks whilst experimenters set up new file on the EEG machine
        end
        PsychPortAudio('Start', tc);
        WaitSecs(0.1);
        exitcode = 0;
        while exitcode==0
            
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, width, height,0)
            textstring = strcat('Click to start block', num2str(currentblock));
            Screen('DrawText', w, textstring, (width/2)-200, (height/2)-400, [0 0 0]);
            
             textstring = strcat('Subject ID number: ', subj);
            Screen('DrawText', w, textstring, (width/2)-300, (height/2)+400, [0 0 0]);               
            
            lastflip = Screen('Flip', w);
            
            [x,y,buttons] = GetMouse;
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if sum(buttons)>0 || keyIsDown
                exitcode = 1;
            end
        end
        
        for n = 1:3                         % 3 trigger pulses to indicate start of block
            if useVP
                Datapixx('SetDoutValues', transformindex(1));
                Datapixx('RegWrRd');
            end
            for nwaitframepairs = 1:6
                Screen('FillRect', w, ST.greylevel);
                drawfixation(w, width, height,0)
                lastflip = Screen('Flip', w);
                
            end
            if useVP
                Datapixx('SetDoutValues', 0);
                Datapixx('RegWrRd');
            end
            for nwaitframepairs = 1:6
                Screen('FillRect', w, ST.greylevel);
                drawfixation(w, width, height,0)
                lastflip = Screen('Flip', w);
                
            end
            
        end
        Screen('FillRect', w, ST.greylevel);
        drawfixation(w, width, height,0)
        lastflip = Screen('Flip', w);
        trialoffset = lastflip;     % last trial was an infinitely long time ago
        trial = 0;
        
        DP.pass = 1;        % always go back to the first pass during the first half of a block
        dpstartingt = sum(squeeze(DP.tno(:,DP.pass)));
        dptrial = dpstartingt;
        
        while trial < E.ntrialsperblock
            trial = trial + 1;
            E.cond = R.triallistordered(1,R.blocktrialsordered(R.trialorder(R.currentrep,trial)));
            E.level = R.triallistordered(2,R.blocktrialsordered(R.trialorder(R.currentrep,trial)));
            R.tno(E.cond,E.level) = R.tno(E.cond,E.level) + 1;
            stimorient = R.orientations(E.cond,E.level,R.tno(E.cond,E.level));
            
            if trial == E.ntrialsperblock/2
                DP.pass = 2;        % change to pass 2 half way through block
            end
            
            targetstim = mkgrating(ST.gratingsize, ST.ncycles, stimorient, E.targetphase, 1) .* window;
            maskstim = mkgrating(ST.gratingsize, ST.maskcycles, stimorient+E.maskori, E.maskphase, 1) .* window;
            
            for n = 1:targetframes              % target and mask
                comp = targetwaveform(n).*targetstim.*E.targetlevelsC(E.level) + maskwaveform(n).*maskstim.*E.masklevelsC(E.level);
                comp = (1+comp)/2;
                if useVP
                    comp = doimagegamma(comp);
                end
                stimlist(n) = Screen('MakeTexture', w, comp, [], [], 2);
            end
            
            nframes = length(stimlist);
            tempdptrial = dptrial;
            temptno = DP.tno;
            for m = 1:DP.lookaheadtrials        % create stimuli for next three trials
                tempdptrial = tempdptrial + 1;
                currenttrial = DP.trialorder(tempdptrial);
                temptno(currenttrial,DP.pass) = temptno(currenttrial,DP.pass) + 1;
                pednoise = DP.seeds(currenttrial,temptno(currenttrial,DP.pass),1);
                targetnoise = DP.seeds(currenttrial,temptno(currenttrial,DP.pass),2);
                testinterval = DP.testintervals(currenttrial,temptno(currenttrial,DP.pass),DP.pass);      % test in interval 1 or 2
                
                for n = 1:dpframes
                    comp = dpwaveform(n).*dpstim.*(targetnoise*DP.noiselevelsC+DP.targetlevelsC(currenttrial));
                    comp = (1+comp)/2;
                    if useVP
                        comp = doimagegamma(comp);
                    end
                    dpstimlist(m,testinterval,n) = Screen('MakeTexture', w, comp, [], [], 2);
                end
                
                for n = 1:dpframes
                    comp = dpwaveform(n).*dpstim.*(pednoise*DP.noiselevelsC);
                    comp = (1+comp)/2;
                    if useVP
                        comp = doimagegamma(comp);
                    end
                    dpstimlist(m,3-testinterval,n) = Screen('MakeTexture', w, comp, [], [], 2);
                end
                
            end
            
            dpcurrentphase = 1;     % 1 is no stimulus, 2 is trial sequence, 3 is wait for response
            dpblankframecount = 1;
            dplookaheadtrial = 0;
            dpwaitperiod = 0;
            Screen('FillRect', w, ST.greylevel);
            drawfixation(w, width, height,1);
            lastflip = Screen('Flip', w, trialoffset+ST.ITI);
            starttime = lastflip;
            
            for n = 1:(ST.nframes*ST.duration)          % loop of one second times ST.duration frames
                
                frameindex = mod(n-1,nframes)+1;
                
                Screen('FillRect',w, ST.greylevel);
                Screen('DrawTextures', w, stimlist(frameindex), [], destRects);
                if dpcurrentphase==2
                    if dpinterval
                        dpframecount = dpframecount + 1;
                        Screen('DrawTexture',w,dpstimlist(dplookaheadtrial,dpinterval,dpframecount), [], dpRect);
                        if dpframecount==1
                            PsychPortAudio('Start', tc);
                        end
                    end
                    if dpframecount == length(dpstimlist)
                        if dpinterval==1
                            dpinterval = 0;
                            dpframecount = 0;
                            dpwaitperiod = 1;
                            dpwaitframes = 0;
                        end
                        if dpinterval==2
                            dpinterval = 0;
                            dpwaitperiod = 0;
                            dpcurrentphase = 3;
                        end
                    end
                    if dpwaitperiod
                        dpwaitframes = dpwaitframes + 1;
                        if dpwaitframes==DP.nframesinISI
                            dpwaitperiod = 0;
                            dpinterval = 2;
                        end
                    end
                end
                drawfixation(w, width, height,1)
                lastflip = Screen('Flip', w, lastflip+ifi*0.5);
                frametime(n) = lastflip - starttime;
                
                if useVP
                    if frameindex==1             % a pulse every 1 second (120 frames)
                        if n==1
                            Datapixx('SetDoutValues', transformindex(E.level*10));     % first trigger also contains condition code
                            Datapixx('RegWrRd');
                        else
                            Datapixx('SetDoutValues', transformindex(1));
                            Datapixx('RegWrRd');
                        end
                    end
                    if frameindex==6           % now set to zero after 3 frame pairs (50ms)
                        Datapixx('SetDoutValues', 0);
                        Datapixx('RegWrRd');
                    end
                end
                
                if dpcurrentphase==3
                    resp = 0;
                    [x,y,buttons] = GetMouse;
                    [keyIsDown, secs, keyCode] = KbCheck;
                    
                    if buttons(1) || keyCode(KbName('LeftArrow'))
                        resp = 1;
                    elseif buttons(2) || keyCode(KbName('RightArrow'))
                        resp = 2;
                    end
                    
                    if resp==testinterval
                        DP.resps(currenttrial,DP.tno(currenttrial,DP.pass),DP.pass) = 1;
                    end
                    if resp>0
                        dpcurrentphase = 1;
                        dpblankframecount = 0;
                    end
                    
                end
                
                if dpcurrentphase==1
                    dpblankframecount = dpblankframecount + 1;
                    if dpblankframecount>=DP.nframesinITI
                        if (ST.duration-(n/ST.nframes))>1.5         % so only do another trial if we have time (1.5 secs)
                            dpcurrentphase = 2;
                            dplookaheadtrial = dplookaheadtrial + 1;
                            dpinterval = 1;
                            dpframecount = 0;
                            dptrial = dptrial + 1;
                            currenttrial = DP.trialorder(dptrial);
                            DP.tno(currenttrial,DP.pass) = DP.tno(currenttrial,DP.pass) + 1;
                            testinterval = DP.testintervals(currenttrial,DP.tno(currenttrial,DP.pass),DP.pass);      % test in interval 1 or 2
                        end
                    end
                end
                
            end
            
            allframetimes(trial,:) = frametime;
            trialoffset = lastflip;
            
            if useVP
                Datapixx('SetDoutValues', transformindex(1));
                Datapixx('RegWrRd');
                for nwaitframepairs = 1:6
                    Screen('FillRect', w, ST.greylevel);
                    drawfixation(w, width, height,0)
                    lastflip = Screen('Flip', w);
                    
                end
                
                Datapixx('SetDoutValues', 0);
                Datapixx('RegWrRd');
            else
                Screen('FillRect', w, ST.greylevel);
                drawfixation(w, width, height,0)
                lastflip = Screen('Flip', w);
            end
            
            if dpcurrentphase==3        % check we're not waiting for a response
                exitcode = 0;
                while ~exitcode
                    resp = 0;
                    [x,y,buttons] = GetMouse;
                    [keyIsDown, secs, keyCode] = KbCheck;
                    
                    if buttons(1) || keyCode(KbName('LeftArrow'))
                        resp = 1;
                        exitcode = 1;
                    elseif buttons(2) || keyCode(KbName('RightArrow'))
                        resp = 2;
                        exitcode = 1;
                    end
                    
                    if resp==testinterval
                        DP.resps(currenttrial,DP.tno(currenttrial,DP.pass),DP.pass) = 1;
                    end
                    
                end
            end
            
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if keyCode(KbName('Escape'))
                breakcode = 1;
                trial = 1000;
                E.block = 100;
                currentblock = 100;
            end
            
            Screen('Close',stimlist(:));
            Screen('Close',dpstimlist(:));
        end
        
        if ~breakcode
            save(fname, 'R', 'E', 'DP');       % save every block of trials
            save(strcat(E.exptpath, 'Results/', E.subj,'TimesBlockM',num2str(R.currentrep),'.mat'), 'allframetimes', 'alltriggertimes');
            load(strcat(E.exptpath,'Results/subjectlist.mat'));
            allblocksdone(sno) = allblocksdone(sno) + 1;
            save(strcat(E.exptpath,'Results/subjectlist.mat'), 'allblocksdone','allIDnos');
        end
        PsychPortAudio('Start', tr);
    end
    
catch
    
    lasterr
    
end

Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

Screen('Flip', w);
Screen('Close',nulltexture);
ShowCursor;

if useVP        % using a ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
end
Screen('CloseAll');
if useVP
    Datapixx('Close');
end

PsychPortAudio('Close', tc);
PsychPortAudio('Close', tr);

end
%--------------------------------------------------------------------------------------------------
function drawfixation(w, width, height, drawframe)

ulx = width/2;
uly = height/2;

Screen('DrawLine', w, [0], ulx-4, uly, ulx+4, uly, 2);
Screen('DrawLine', w, [0], ulx, uly-4, ulx, uly+4, 2);

if drawframe
    r = [0 0 110 110];
    r = CenterRectOnPoint(r, ulx, uly);
    Screen('FrameRect', w, [0], r, 2);
end

% if drawframe
%     imsize = 110;
%     r = [0 0 imsize+2 imsize+2];
%     for x = 1:9
%         for y = 1:5
%             offsetx = imsize*(x-5);
%             offsety = imsize*(y-3);
%             r = CenterRectOnPoint(r, ulx+offsetx, uly+offsety);
%             Screen('FrameRect', w, [0], r, 2);
%         end
%     end
% end

end
%--------------------------------------------------------------------------
function mouseloop

exitcode = 0;

while exitcode==0
    [x,y,buttons] = GetMouse;
    
    if sum(buttons)>0
        exitcode = 1;
    end
end

end
%--------------------------------------------------------------------------
function dB = d_PercentTodB(perc)

%converts from % contrast to contrast in dB

dB = 20 * log10(perc);

end
%--------------------------------------------------------------------------------------------------
function [perc] = d_dBtoPercent(dB)

%converts from dB to % contrast

perc = 10.^(dB/20);

end
%--------------------------------------------------------------------------
function kbloop

exitcode = 0;

while exitcode==0
    [keyIsDown, secs, keyCode] = KbCheck;       % also monitor keyboard for breaks
    
    if sum(keyCode(79:80))>0        % respond to left and right arrows
        exitcode = 1;
    end
end

end
%--------------------------------------------------------------------------------------------------
function output = doimagegamma(i)

% gamma corrects the stimuli before sending to Bits++
% adapted from Mark's code, DHB 29.01.08

% parameters from last gamma correct, 12/8/14 on VPixx in M16 mode
k = 0.5344;
Lmax = 101.8721;
j0 = 0.1292;
gamma = 1.9252;
%%%%%

i0 = 0;
imax = 1;                                               % Bits++ always scaled between 0 and 1
imean = (i0+imax)/2;
jmax = 1;

Lmin = k + (Lmax-k)*(max(-j0,0)/(jmax-j0) ).^gamma;     % Eqn 2, with j set to 0, to get Lmin
Lmin = max(Lmin,0);                                     % ensure Lmin not <0
Lmean = (Lmin+Lmax)/2;
L = Lmean + (Lmax-Lmean)*(i-imean)/(imax-imean);        % desired luminance values Eqn 4
j = ((L - k)/(Lmax-k)).^(1/gamma)*(jmax - j0) + j0;     % These are the gamma-corrected lut values, j: Eqn 3
output = max(j,j0);                                     % Eqn 3 conditional
output = double(output);

end
%--------------------------------------------------------------------------
function imag1 = mkgrating(Regionsize, f, o, p, c)

%TSM; 26.6.03
% modified by DHB to make single component gratings only, scaled from -1 to 1
% f is spatial frequency, scaled as cycles per image
% o is orientation (degrees), p is phase (degrees relative to centre), c is contrast

p = p*pi/180;
o = o*2*pi/360;		% convert from degrees to radians
f = f/Regionsize;
x0 = ((Regionsize+1)/2);
y0 = x0;

u = f .* cos(o) * 2 * pi;
v = f .* sin(o) * 2 * pi;

imag1 = zeros(Regionsize, Regionsize);
[xx, yy] = meshgrid(1:Regionsize, 1:Regionsize);

imag1(:,:) = (c .* sin(u .*(xx-x0) + v.*(yy-y0) + p));

end
%--------------------------------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% ends an array 'mask' that is 1 inside the circular window, shading to zero outside
% W, H are the width and height of the whole (rectangular or square) array, in pixels
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is the smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies this diameter (in relative units, range 0 -> 1)
% MAG, 27.2.04

%soft window parameters
if nargin<3, D = 0.9; end % sets default diameter to 0.9
radius = min(W*D/2,H*D/2);% radius in pixels
blur = 2*(min(W/2,H/2) - radius);  % blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
%image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% make circular soft window
mask = single((xx.*xx + yy.*yy) < radius^2); % logical 0 outside the circle,1 inside it
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask/max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%--------------------------------------------------------------------------
function output = transformindex(input)

% fixes the binary inputs for the EEG amplifier because the pins are in a different order from the ViewPixx
% desired numbers must be <256
% DHB 18/8/14

truebits = 2.^(2:2:24);
dn = dec2bin(input,length(truebits));
output = 0;
for m = 1:length(truebits)
    output = output + truebits(m)*str2num(dn(end-m+1));
end

end
%--------------------------------------------------------------------------
function [subj, blocksdone, sno] = runsubjectroutine(exptpath)

load(strcat(exptpath,'Results/subjectlist.mat'));

textstring = 'Please enter subject ID number';

exitcode = 0;
while ~exitcode
IDnumber = inputdlg(textstring);
IDnumber = str2num(IDnumber{1});
a = find(allIDnos==IDnumber);

if ~isempty(a)      % subject already exists
   sno = a;
   blocksdone = allblocksdone(a);
else
    sno = length(allIDnos)+1;
    allIDnos(end+1) = IDnumber;
    blocksdone = 0;
end
    
if blocksdone==4        % someone already ran this
    textstring = 'Number in use. Please enter an ID number';
else
    exitcode = 1;
end
end

subj = strcat('S',num2str(IDnumber));
% if length(subj)<2
%     subj = strcat('0',subj);
% end
% if length(subj)<3
%     subj = strcat('0',subj);
% end
% subj = strcat('S',subj);

save(strcat(exptpath,'Results/subjectlist.mat'), 'allblocksdone','allIDnos');

end
%--------------------------------------------------------------------------