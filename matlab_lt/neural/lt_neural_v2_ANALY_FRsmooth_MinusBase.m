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
% changed to be before onset of WN...


%%
%% ####################### FOR EACH MOTIF/NEURON, SUBTRACT BASELINE SM FR
lt_figure; hold on;
title('baseline and last epoch data windows');
ylaball = {};

numbirds = length(SwitchStruct.bird);
maxneur = max(OUTDAT.All_neurnum);
ycount = 1;

AllBase_indsepoch = cell(size(OUTDAT.All_FRsmooth,1),1);
AllWN_indsepoch = cell(size(OUTDAT.All_FRsmooth,1),1);

AllMinusBase_FRmeanAll = cell(size(OUTDAT.All_FRsmooth,1),1);
AllMinusBase_FRsemAll = cell(size(OUTDAT.All_FRsmooth,1),1);
AllMinusBase_tbinAll = cell(size(OUTDAT.All_FRsmooth,1),1);
AllBase_FRsmooth_Hi = cell(size(OUTDAT.All_FRsmooth,1),1);
AllBase_FRsmooth_Lo = cell(size(OUTDAT.All_FRsmooth,1),1);
AllMinusBase_PitchZ = nan(size(OUTDAT.All_FRsmooth,1), length(prctile_divs));

AllTargLearnDir = nan(size(OUTDAT.All_FRsmooth,1),1);

for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            indstmp = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss);
            if ~any(indstmp)
                continue
            end
            
            swtime = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            ldir = cell2mat(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(2:2:end));
            ldir = unique(ldir);
            assert(length(ldir)==1);
            
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if ~any(indsthis)
                    continue
                end
                
%                 assert(length(indsthis)==1);
                
                % ================ COLLECT FRATE
                for j=indsthis'
                    
                    if ss==1 & nbasetime_ignoreswitch1==1
                        % then ignore base time
                        [FRmeanAll, FRsemAll, tbin, frlo, frhi, PitchmeanAll, TimesAll, ...
                            timesbase, ffbase, indsbase_epoch, indsWN_epochOnsets] = ...
                            fn_subtractbase(OUTDAT, j, prctile_divs, usepercent);
                        
                    else
                        [FRmeanAll, FRsemAll, tbin, frlo, frhi, PitchmeanAll, TimesAll, ...
                            timesbase, ffbase, indsbase_epoch, indsWN_epochOnsets] = ...
                            fn_subtractbase(OUTDAT, j, prctile_divs, usepercent, ...
                            nbasetime, swtime);
                    end
                    
                    
                    if (0) % sanity check, make sure inds are what I think they are.
                        tmp = mean(OUTDAT.All_FRsmooth{j,2}(:, indsWN_epochOnsets(end):end),2) - ...
                            mean(OUTDAT.All_FRsmooth{j,1}(:,indsbase_epoch),2)
                        %                     tmp = cellfun(@(x)mean(x,2), OUTDAT.All_FRsmooth(j,:), 'UniformOutput', 0);
                        figure; hold on;
                        plot(tmp, 'k');
                        %                     plot(tmp{2}-tmp{1}, 'k');
                        plot(FRmeanAll{3}, 'r');
                    end
                    
                    % =========== SAVE
                    assert(isempty(AllBase_indsepoch{j}));
                    
                    AllTargLearnDir(j) = ldir;
                    AllBase_indsepoch{j} = indsbase_epoch;
                    AllWN_indsepoch{j} = indsWN_epochOnsets;
                    AllMinusBase_FRmeanAll{j} = FRmeanAll;
                    AllMinusBase_FRsemAll{j} = FRsemAll;
                    AllMinusBase_tbinAll{j} = tbin;
                    AllBase_FRsmooth_Hi{j} = frhi;
                    AllBase_FRsmooth_Lo{j} = frlo;
                    AllMinusBase_PitchZ(j, :) = PitchmeanAll;
                end
            end
            
            % ========== COLLECT TO PLOT ACROSS SWITCHES
            disp('NOTE: This might not be perfectly accurate if different neurons for this switch have different duration data');
            disp('In that case will only thne reflect timing for the last neuron iterated over');
            
            timetolock = timesbase(end);
            tbasetmp = (timesbase-timetolock)*24;
            line([min(tbasetmp) max(tbasetmp)], [ycount ycount], 'LineWidth', 2, 'Color','r');
            
            twntmp = (TimesAll{end}-timetolock)*24;
            %             ffzwntmp = PitchmeanAll(end)*learndir;
            line([min(twntmp) max(twntmp)], [ycount ycount], 'LineWidth', 2, 'Color', 'r');
            
            line([max(tbasetmp) min(twntmp)], [ycount ycount], 'Color', [0.7 0.7 0.7], 'LineStyle', '--')
            
            
            % --- note down expt
            strthis = [birdname '-' exptname(end-7:end) '-sw' num2str(ss)];
            %             lt_plot_text(max(twntmp)+0.1, ycount, strthis, 'k');
            ylaball = [ylaball; strthis];
            ycount = ycount+1;
            
        end
        line([xlim], [ycount-0.5 ycount-0.5], 'Color', [0.3 0.3 0.7]);
    end
end

OUTDAT.AllMinusBase_FRmeanAll = AllMinusBase_FRmeanAll;

OUTDAT.AllBase_indsepoch = AllBase_indsepoch;
OUTDAT.AllWN_indsepoch= AllWN_indsepoch;

OUTDAT.AllMinusBase_FRsemAll = AllMinusBase_FRsemAll;
OUTDAT.AllMinusBase_tbinAll = AllMinusBase_tbinAll;
OUTDAT.AllBase_FRsmooth_Hi = AllBase_FRsmooth_Hi;
OUTDAT.AllBase_FRsmooth_Lo = AllBase_FRsmooth_Lo;
OUTDAT.AllMinusBase_PitchZ = AllMinusBase_PitchZ;

OUTDAT.AllTargLearnDir = AllTargLearnDir;

set(gca, 'YTick', 1:ycount-1);
set(gca, 'YTickLabel', ylaball);
end



function [FRmeanAll, FRsemAll, tbin, frlo, frhi, PitchmeanAll, TimesAll, ...
    timesbase, ffbase, indsbase_epoch, indsWN_epochOnsets] = fn_subtractbase(OUTDAT, j, ...
    prctile_divs,usepercent, nbasetime, swtime)
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
    swtime = 0;
end

% ====== RUN
% =========== 1) collect baseline mean (use last N trials)
indtmp = 1;

% --- minimum time for base trials
ttmp = OUTDAT.All_FF_t{j,2}(1);
% mintime = ttmp - nbasetime/(24*60); % the earliest data to take for baseline
mintime = swtime - nbasetime/(24*60);

% -- which bnase inds to keep?
tbase = OUTDAT.All_FF_t{j,1};
indstokeep = tbase>=mintime;
indsbase_epoch = find(indstokeep);

% -- get mean FR
frmat = OUTDAT.All_FRsmooth{j, indtmp}(:, indstokeep);
ff = OUTDAT.All_FF{j, indtmp}(indstokeep);

frmean_base = mean(frmat, 2);
frbase_std = std(frmat, [], 2);
ff_base = mean(ff);
ff_base_std = std(ff);

% -------------- OUTPUT TIME OF ALL BASE INDS
timesbase = tbase(indstokeep);
ffbase = ff_base;

if (0)
    
    tbase = OUTDAT.All_FF_t{j, 1};
    twn = OUTDAT.All_FF_t{j,2};
    fbase = OUTDAT.All_FF{j,1};
    fwn = OUTDAT.All_FF{j,2};
    figure; hold on;
    
    if isnan(fbase(1))
        plot(tbase, 1, 'ok');
        plot(twn, 2, 'or');
    else
        plot(tbase, fbase, 'ok');
        plot(twn, fwn, 'or');
    end
end
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
inddivs(end) = inddivs(end)+1; % since want to end on last trial.
PitchmeanAll = nan(1, length(prctile_divs));

TimesAll = cell(1, length(prctile_divs));
for k=1:length(inddivs)-1
    
    indsthis = inddivs(k):inddivs(k+1)-1;
    
    % ========= FRATE
    ymean = mean(frmat(:, indsthis), 2);
    ystd = std(frmat(:, indsthis), [], 2);
    ysem = lt_sem(frmat(:, indsthis)')';
    
    if usepercent==1
        ymean = 100*(ymean - frmean_base)./frmean_base;
    elseif usepercent==0
        ymean = ymean - frmean_base;
    elseif usepercent==2
        % then zscore
        ymean = (ymean-frmean_base)./frbase_std;
    elseif usepercent==3
        ymean = (ymean-frmean_base)./sqrt(0.5.*((frbase_std.^2)+(ystd.^2)));
    end
    
    % =========== FF
    ffmean = mean(ffmat(indsthis));
    ffmean = (ffmean - ff_base)./ff_base_std;
    
    TimesAll{k} = tvals(indsthis); % collects times for last bin.
    
    % ========= OUTPUT
    if any(isnan(ymean))
        keyboard
    end
    
    FRmeanAll{k} = ymean;
    FRsemAll{k} = ysem;
    PitchmeanAll(k) = ffmean;
end

% --------- OUTPUT TIME OF TRAINING EPOCH INDS
indsWN_epochOnsets = inddivs;

end


