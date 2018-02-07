function testtriggers

    Datapixx('Open');
    Datapixx('StopAllSchedules');
    Datapixx('RegWrRd');    % Synchronize Datapixx registers to local register cache
    
    % We'll make sure that all the TTL digital outputs are low before we start
    Datapixx('SetDoutValues', 0);
    Datapixx('RegWrRd');
    n = 0;
    while n < 10
        n = n + 1
            Datapixx('SetDoutValues', transformindex(n));
            Datapixx('RegWrRd');
            WaitSecs(0.2);
            val = Datapixx('GetDoutValues')
%             breakcode = 0;
%             while ~breakcode
%             [keyIsDown, secs, keyCode] = KbCheck;
%             if keyCode(KbName('RightArrow'))
%                 breakcode = 1;
%             end
%             if keyCode(KbName('Escape'))
%                 n = 1000;
%                 breakcode = 1;
%             end
%             end
    end
            
            Datapixx('Close');
            
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