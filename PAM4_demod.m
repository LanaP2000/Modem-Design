% This function will do 4-PAM demodulation
% It's will demapping for non-Gray mapping
%   -3 into 00 
%   -1 into 01
%    1 into 10
%    3 into 11

function bits = PAM4_demod(code)
bits = zeros(1,2*length(code));

% the following code is to do the demodulation, the code is written for
% easy reproduction in your project. Efficient coding can be done with 1 or
% 2 MATLAB lines.

for i=1:length(code)
    if code(i) < -2
        bits(2*i-1:2*i) = [0 0];
    elseif code(i) < 0
        bits(2*i-1:2*i) = [0 1];
    elseif code(i) < 2
        bits(2*i-1:2*i) = [1 0];
    else
        bits(2*i-1:2*i) = [1 1];
    end
end
