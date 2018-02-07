function output = transformindex(input)

truebits = 2.^(2:2:24); % 4, 16, 64, 256, 1024, 4096, 16384... 16777216. Basically 2 ^ 2/4/6/8...24. 

dn = dec2bin(input,length(truebits)); % transforms input to a binary representation as a string with at least '12' bits. So 10 becomes 00000000001010. 

output = 0;

for m = 1:length(truebits)
    
    output = output + truebits(m) * str2num(dn(end-m+1));
    
end

end

% for i = 1:10
% transformindex(i) = 4, 16, 20, 64, 68, 80, 84, 256, 260, 272. 
% end