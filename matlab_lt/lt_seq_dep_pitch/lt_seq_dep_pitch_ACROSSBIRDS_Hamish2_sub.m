function [ffmean, ffcv, ffvals] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
    birdnum, exptnum, daystoget, syltoget, doMUSC)
%% lt 7/9/18 - extracts raw FF during PBS and MUSC, mean of day means

i = birdnum;
ii = exptnum;

%%
DatThis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning;

if doMUSC==0
% ============ PBS
    t_byday = DatThis.AllDays_PlotLearning.DataMatrix.(syltoget).Tvals_WithinTimeWindow(daystoget);
    ff_byday = DatThis.AllDays_PlotLearning.DataMatrix.(syltoget).FFvals_WithinTimeWindow(daystoget);

elseif doMUSC==1
   % =========== MUSC 
    t_byday = DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(syltoget).Tvals_WithinTimeWindow(daystoget);
    ff_byday = DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(syltoget).FFvals_WithinTimeWindow(daystoget);
end

% ============ extract mean of day stats
    ffmean = mean(cellfun(@mean, ff_byday)); % mean
    ffcv = mean(cellfun(@std, ff_byday)./cellfun(@mean, ff_byday)); % cv
    ffvals = cell2mat(ff_byday); % ffvals

