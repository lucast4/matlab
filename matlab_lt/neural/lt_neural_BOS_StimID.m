function [datdig, idthis] = lt_neural_BOS_StimID(datdig, fs)
%% lt 11/30/18 - given digital signal (syls) appends marker so can identify BOS type.

% -- places 1ms (2ms period) blips starting 0.1 sec from song offset.
idthis = randi(100,1);

% --- place signal into dig signal
% find the last offset, then place this signal one second after
[~, ftimes] = falltime(datdig);

% - 1ms on, 1ms off
onstmp = 0:(0.002*fs):(0.002*fs)*(idthis-1);
offstmp = onstmp+0.001*fs;
% -- add them to the last offset + 0.1 sec
onstmp = ceil(ftimes(end))+0.1*fs+onstmp;
offstmp = ceil(ftimes(end))+0.1*fs+offstmp;

% --- make sure there is enough 0 data
assert(length(datdig)>offstmp(end))
% --- add them to signal
for j=1:length(onstmp)
    datdig(onstmp(j):offstmp(j)) = 1;
end
