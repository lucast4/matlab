function [Stats] = lt_neural_BOS_StatsUnits(AllUnits, SummaryBOS, PARAMS)
%% lt 12/3/18 - Gets stats across units (e.g. mean resposne
% input AllUnits must already be filtered. will assumes want to get stats
% across all units in AllUnits.

% ============= CHECK THAT EACH UNIT ONLY OCCUPIES ONE CELL.
units = AllUnits.unitnum;
expt = AllUnits.exptnum;
units_unique = lt_tools_grp2idx({expt, units});
assert(max(units_unique)==length(units_unique), 'some unts have multipel cells.., means did not filter correctly the inputs to this function');

%% GET MEAN RESPONSE (RATE VS.TIME)
% ---- for each unit, get mean response to each mtoif
frmat = AllUnits.FRmat;
frmat_t = AllUnits.FRmat_t;

% --- GET MEAN RESPONSE FOR EACH UNIT/MOTIF
frmean = cellfun(@(x)mean(x,2), frmat, 'UniformOutput', 0);
frsem = cellfun(@(x)lt_sem(x'), frmat,  'UniformOutput', 0);
tmean = cellfun(@(x)mean(x,2), frmat_t, 'UniformOutput', 0);

% --- TAKE MEAN ACROSS UNITS
ymean = mean(cell2mat(cellfun(@transpose, frmean, 'UniformOutput', 0)));
ysem = lt_sem(cell2mat(cellfun(@transpose, frmean, 'UniformOutput', 0)));
t = mean(cell2mat(cellfun(@transpose, tmean, 'UniformOutput', 0)));

% =================== OUTPUT
Stats.frate_unitmean = frmean;
Stats.frate_unitsem = frsem;
Stats.frate_mean = ymean;
Stats.frate_sem = ysem;
Stats.frate_t = t;


%% GET MEDIAN SYL ONSET, OFFSET
on = AllUnits.SylOnsets;
off = AllUnits.SylOffsets;

% --- GET MEAN RESPONSE FOR EACH UNIT/MOTIF
on = cellfun(@(x)median(x,1), on, 'UniformOutput', 0); % across trials
on = median(cell2mat(on)); % across units.

off = cellfun(@(x)median(x,1), off, 'UniformOutput', 0); % across trials
off = median(cell2mat(off)); % across units.

Stats.sylon_median = on;
Stats.syloff_median = off;
