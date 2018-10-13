function [SwitchTimesMat, SwitchTimesMat_conting] = lt_tools_SwitchTimes(SwitchTimes, lastdatenum)
%% lt 10/11/18 - given SwitchTimes cell, outputs mat with epochs

% e.g. cell which has date (in ddmmmyyyy-HHMM) format, followed by 2-length
% strings with conditions (i.e. up, of, down, for WN - currently these are
% the only ones allowed)

%  INPUTS
% SwitchTimes = {...
%     '29Sep2018-2021', 'of-up', ...
%     '04Oct2018-2400', 'up-dn', ...
%     '09Oct2018-2400', 'dn-up'...
%     };

% lastdatenum, datenum to terminate the last epoch on.

% OUTPUT, matrices with size in dim 1 one larger than inpu tbeucase is for
% p[hases, not switches. 

% SwitchTimesMat =
% 
%                          0          737332.847916667
%           737332.847916667                    737338
%                     737338                    737343
%                     737343                    737343
% 
% SwitchTimesMat_conting =
% 
%      0
%      1
%     -1
%      1

% ------ convert switchtimes to range of datetimes
SwitchTimesMat = nan(length(SwitchTimes)/2+1, 2);
SwitchTimesMat_conting = nan(length(SwitchTimes)/2+1,1);
for j=1:length(SwitchTimes)/2+1
   if j==1
       % -- begin to now
       t1 = 0;
       t2 = datenum(SwitchTimes{(j)*2-1}, 'ddmmmyyyy-HHMM'); % current switch
       conting = SwitchTimes{(j)*2}(1:2);
   elseif j==length(SwitchTimes)/2+1
       % -- to end
       t1 = datenum(SwitchTimes{(j-1)*2-1}, 'ddmmmyyyy-HHMM'); % previso switch
       t2 = lastdatenum;
       conting = SwitchTimes{(j-1)*2}(4:5);
   else
       t1 = datenum(SwitchTimes{(j-1)*2-1}, 'ddmmmyyyy-HHMM'); % previso switch
       t2 = datenum(SwitchTimes{(j)*2-1}, 'ddmmmyyyy-HHMM'); % current switch
       conting = SwitchTimes{(j-1)*2}(4:5);
   end
  
   SwitchTimesMat(j, :) = [t1 t2];
   if strcmp(conting, 'of')
       c = 0;
   elseif strcmp(conting, 'dn')
       c = -1';
   elseif strcmp(conting, 'up')
       c = 1;
   end
   SwitchTimesMat_conting(j) = c;
end
