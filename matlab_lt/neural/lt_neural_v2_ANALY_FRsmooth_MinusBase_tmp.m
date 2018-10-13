function OUTDAT = lt_neural_v2_ANALY_FRsmooth_MinusBase(OUTDAT, SwitchStruct,...
    usepercent, nbasetime, nbasetime_ignoreswitch1, prctile_divs)
%% lt 9/11/18 - subtracts baseline smoothed FR from fr during training.

% prctile_divs = [33 66 100]; % percentiles to divide up data by
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
% prctile_divs = [50 100]; % percentiles to divide up data by
% epochtoplot = 2; % i.e. out of the epochs decided by prctile_divs

% usepercent = 0; % if 1, then gets fr percent deviation from baseline. otherwise uses hz diff
% nbasetime = 60; % 60 minutes bnefore first WN trial is min time for baseline
% nbasetime = []; % 60 minutes bnefore first WN trial is min time for baseline


%%

%% ####################### FOR EACH MOTIF/NEURON, SUBTRACT BASELINE SM FR
lt_figure; hold on;


numbirds = length(SwitchStruct.bird);
maxneur = max(OUTDAT.All_neurnum);

AllMinusBase_FRmeanAll = cell(size(OUTDAT.All_FRsmooth,1),1);
AllMinusBase_FRsemAll = cell(size(OUTDAT.All_FRsmooth,1),1);
AllMinusBase_tbinAll = cell(size(OUTDAT.All_FRsmooth,1),1);
AllBase_FRsmooth_Hi = cell(size(OUTDAT.All_FRsmooth,1),1);
AllBase_FRsmooth_Lo = cell(size(OUTDAT.All_FRsmooth,1),1);
AllMinusBase_PitchZ = nan(size(OUTDAT.All_FRsmooth,1), length(prctile_divs));
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if ~any(indsthis)
                    continue
                end
                
                % ================ COLLECT FRATE
                for j=indsthis'
                    
                    if ss==1 & nbasetime_ignoreswitch1==1
                        % then ignore base time
                        [FRmeanAll, FRsemAll, tbin, frlo, frhi, PitchmeanAll, TimesAll, ...
                            timesbase, ffbase] = fn_subtractbase(OUTDAT, j, prctile_divs, usepercent, [])
                    else
                        [FRmeanAll, FRsemAll, tbin, frlo, frhi, PitchmeanAll, TimesAll, ...
                            timesbase, ffbase] = fn_subtractbase(OUTDAT, j, prctile_divs, usepercent, nbasetime)
                    end
                    
                    % =========== SAVE
                    AllMinusBase_FRmeanAll{j} = FRmeanAll;
                    AllMinusBase_FRsemAll{j} = FRsemAll;
                    AllMinusBase_tbinAll{j} = tbin;
                    AllBase_FRsmooth_Hi{j} = frhi;
                    AllBase_FRsmooth_Lo{j} = frlo;
                    AllMinusBase_PitchZ(j, :) = PitchmeanAll;
                end
            end
            
            % ========== COLLECT TO PLOT ACROSS SWITCHES
            
        end
    end
end

OUTDAT.AllMinusBase_FRmeanAll = AllMinusBase_FRmeanAll;
OUTDAT.AllMinusBase_FRsemAll = AllMinusBase_FRsemAll;
OUTDAT.AllMinusBase_tbinAll = AllMinusBase_tbinAll;
OUTDAT.AllBase_FRsmooth_Hi = AllBase_FRsmooth_Hi;
OUTDAT.AllBase_FRsmooth_Lo = AllBase_FRsmooth_Lo;
OUTDAT.AllMinusBase_PitchZ = AllMinusBase_PitchZ;

end



function [FRmeanAll, FRsemAll, tbin, frlo, frhi, PitchmeanAll, TimesAll, ...
    timesbase, ffbase] = fn_subtractbase(OUTDAT, j, prctile_divs, usepercent, nbasetime)
% ===== PARAMS
% j is ind to extract from in OUTDAT.
% nbasetime in minutes, counting back from first WN song

% --- OUT:
% FRmeanAll, each cell is one period during learning, mean smoothed fr
% minus baseline.

if ~exist('nbasetime', 'var')
    nbasetime = [];
end
if isempty(nbasetime)
    nbasetime = 10000000;
end
% ====== RUN
% =========== 1) collect baseline mean (use last N trials)
indtmp = 1;

% --- minimum time for base trials
ttmp = OUTDAT.All_FF_t{j,2}(1);
mintime = ttmp - nbasetime/(24*60); % the earliest data to take for baseline

% -- which bnase inds to keep?
tbase = OUTDAT.All_FF_t{j,1};
indstokeep = tbase>=mintime;

% -- get mean FR
frmat = OUTDAT.All_FRsmooth{j, indtmp}(:, indstokeep);
ff = OUTDAT.All_FF{j, indtmp}(indstokeep);

frmean_base = mean(frmat, 2);
ff_base = mean(ff);
ff_base_std = std(ff);

% -------------- OUTPUT TIME OF ALL BASE INDS
timesbase = tbase(indstokeep);
ffbase = ff_base;

% ------------- collect mean baseline divided by median of FF
if (0)
    thresh = median(ff);
    frlo = mean(frmat(:, ff<thresh), 2);
    frhi = mean(frmat(:, ff>thresh), 2);
else
    frlo = mean(frmat(:, ff<prctile(ff, [33.3])), 2);
    frhi = mean(frmat(:, ff>prctile(ff, [66.6])), 2);
end

% ============= 2) during training, bin by time period
indtmp = 2;

frmat = OUTDAT.All_FRsmooth{j, indtmp};
tbin = OUTDAT.All_FRsmooth_t{j, indtmp};
tvals = OUTDAT.All_FF_t{j, indtmp};
ffmat = OUTDAT.All_FF{j, indtmp};

% ---- bin
FRmeanAll = cell(1, length(prctile_divs));
FRsemAll = cell(1, length(prctile_divs));
inddivs = [1 ceil(prctile(1:size(frmat,2), prctile_divs))]; % get inds that divide up bins.
PitchmeanAll = nan(1, length(prctile_divs));

TimesAll = cell(1, length(prctile_divs));
for k=1:length(inddivs)-1
    
    % ========= FRATE
    ymean = mean(frmat(:, inddivs(k):inddivs(k+1)), 2);
    ysem = lt_sem(frmat(:, inddivs(k):inddivs(k+1))')';
    
    if usepercent==1
        ymean = 100*(ymean - frmean_base)./frmean_base;
    else
        ymean = ymean - frmean_base;
    end
    
    % =========== FF
    ffmean = mean(ffmat(inddivs(k):inddivs(k+1)));
    ffmean = (ffmean - ff_base)./ff_base_std;
    
    TimesAll{k} = tvals(inddivs(k):inddivs(k+1)); % collects times for last bin.
    
    % ========= OUTPUT
    if any(isnan(ymean))
        keyboard
    end
    
    FRmeanAll{k} = ymean;
    FRsemAll{k} = ysem;
    PitchmeanAll(k) = ffmean;
end

% --------- OUTPUT TIME OF TRAINING EPOCH INDS


end


