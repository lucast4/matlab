function [OUTSTRUCT, PARAMS, LFPXCORR_Base_tr_f_cclag, LFPXCORR_WN_tr_f_cclag, ...
    LFPXCORR_freqsall] = lt_neural_LFP_RecalcCoh(OUTSTRUCT, SwitchCohStruct, LFPSTRUCT, SwitchStruct, ...
    PARAMS, usecorrcoeff, thingstodo)
%% lt 12/14/18 - recalcualte things like coherence, lfp xcorr...

assert(length(thingstodo)==1, 'currently the last thing will override all other thiongs ...');

%% ######################## LFP XCORR PARAMS
% Unbiased, but not normal;ized (i.e. is not corr coeff. this 
% NOTE: takes xcorr of the real part of hilvert trasnform.
lfpxcorr_dir = ['/bluejay5/lucas/analyses/neural/FILTER/MOTIFEXTRACT/' PARAMS.savemarker];
lfpxcorr_twindtoget = [-0.1 0.02]; % relative syl onset
lfpxcorr_ffwindtoget = [20 40]; % inclusive [min max]
% maxlag = 0.05*1500; % 50 ms % NOW is deteriend automatically based on
% peridicyt.

usecorrcoeff = 1;

LFPXCORR_Base_tr_f_cclag = cell(size(OUTSTRUCT.bnum,1),1);
LFPXCORR_WN_tr_f_cclag = cell(size(OUTSTRUCT.bnum,1),1);
LFPXCORR_freqsall = cell(size(OUTSTRUCT.bnum,1),1);

fs = 1500;

onlyOnePeriod =1 ; % then only looks at peaks less than one period off from 0 (xcorr lags);
% will make it be 75% of the periopd.

%% ######################## COHERENCE PARAMS
% ==== 1) USE ACROSS TRIALS X TAPERS,
% equalizeN = 1; % trial number. bootstraps muyltipel subsamples.

% cohversion = 'mtaper_all'; % trials x tapers number of datapoints..
% cohversion = 'mtaper_trials'; % one coh for each trial, then takes average
cohversion = 'welch_trials'; % one coh for each trial, then takes average [divides up each segment into 1/2 length windows, and gets 4 overlapping ones]

trialstouse = 'default'; % then inds_epoch as in origianl analysis
% trialstouse = 'default(matched)'; % then inds_epoch as in origianl analysis, but downsamples whichever one is longer.
% trialstouse = 'onlyescapes(matched)'; % then only escapes (entire training duration), sample sizes matched
% trialstouse = 'onlyescapes'; % escapes, takes all trials ignore sample size matching
% trialstouse = 'onlyescapes (epoch)'; % escapes, takes only trials within epoch (ignore sample size matching)
% trialstouse = 'all'; % entire wn vs. base
% (takes entire train and base).
ntapers = []; % leave [] to set as default.
movingwin = [0.1 0.01]; % leave exmpty for default. [applies for both welches and multiutaper]
tw = [];

CohMean_Base = cell(size(OUTSTRUCT.bnum,1),1);
CohMean_WN = cell(size(OUTSTRUCT.bnum,1),1);
Tall = [];
Fall = [];


%% RUN
SampleSizes = []; % base wn
for i=1:length(SwitchCohStruct.bird)
    for ii=1:length(SwitchCohStruct.bird(i).exptnum)
        
        for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
            
            for mm=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                                    disp([i ii ss mm]);

                datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm);
                if isempty(datthis.bregionpair)
                    continue
                end
                pairstoget = datthis.chanpair;
                tvals = datthis.tvals;
                setthis = datthis.neursetused;
                
                t_LFP = LFPSTRUCT.bird(i).experiment(ii).setnum(setthis).motif(mm).t_relons;
                
                swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
                swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
                swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
                
                
                ishit = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).WNhit;
%                 ishit_time = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).WNhittimes_min;
                
                if strcmp(trialstouse, 'default')
                    inds_base = datthis.indsbase_epoch;
                    inds_WN = datthis.indsWN_epoch;
                elseif strcmp(trialstouse, 'default(matched)')
                    inds_base = datthis.indsbase_epoch;
                    inds_WN = datthis.indsWN_epoch;
                    ntokeep = min([length(inds_base) length(inds_WN)]);
                    inds_base = inds_base(end-ntokeep+1:end);
                    inds_WN = inds_WN(end-ntokeep+1:end);
                elseif strcmp(trialstouse, 'onlyescapes(matched)')
                    % --- get wn and base trials that aren't hit
                    inds_WN = find(tvals>swthis & tvals<swpost & ishit==0);
                    inds_base = find(tvals<swthis & tvals>swpre & ishit==0);
                    
                    ntrials = min([length(inds_base) length(inds_WN)]);
                    
                    % --- prune sample size to match
                    inds_base = inds_base(end-ntrials+1:end);
                    inds_WN = inds_WN(end-ntrials+1:end);
                elseif strcmp(trialstouse, 'onlyescapes')
                    % --- get wn and base trials that aren't hit
                    inds_WN = find(tvals>swthis & tvals<swpost & ishit==0);
                    inds_base = find(tvals<swthis & tvals>swpre & ishit==0);
                elseif strcmp(trialstouse, 'all')
                    inds_WN = find(tvals>swthis & tvals<swpost);
                    inds_base = find(tvals<swthis & tvals>swpre);
                elseif strcmp(trialstouse, 'onlyescapes (epoch)')
                    inds_base = intersect(find(ishit==0), datthis.indsbase_epoch);
                    inds_WN = intersect(find(ishit==0), datthis.indsWN_epoch);
                else
                    assert(1==2, 'need to code..');
                end
                
                SampleSizes = [SampleSizes; [length(inds_base) length(inds_WN)]];
                
                % =============== LOAD LFP FILTER BANK, IF DESIRED
                if any(strcmp(thingstodo, 'lfpxcorr'))
                    % ---- load
                    fname = [lfpxcorr_dir '/filtdat_bird' num2str(i) ...
                        '_expt' num2str(ii) '_set' num2str(setthis) '_mot' num2str(mm) '.mat'];
                    filtdatstruct = load(fname);
                    filtdatstruct = filtdatstruct.filtdatstruct;
                end
                
                % ============ FOR EACH PAIR, GET BASELINE AND WN COHERENCE
                npairs = size(pairstoget,1);
                for cc=1:npairs
                    chan1 = pairstoget(cc,1);
                    chan2 = pairstoget(cc,2);
                    
                    lfp1 = datthis.lfpall(:, datthis.lfpall_chans==chan1);
                    lfp2 = datthis.lfpall(:, datthis.lfpall_chans==chan2);
                    
                    
                    % --- find location of this in OUTSTRUCT
                    indthis_out = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
                        OUTSTRUCT.motifnum==mm & OUTSTRUCT.chanpair(:,1)==chan1 ...
                        & OUTSTRUCT.chanpair(:,2)==chan2);
                    
                    if isempty(indthis_out)
                        disp('missing this motif, skip');
                        continue
                    end
                    assert(length(indthis_out)==1, 'not in outstruct..');
                    
                    
                    % ##########################  BASELINE
                    indstoget = inds_base;
                    
                    lfp1_tmp = cell2mat(cellfun(@transpose, lfp1(indstoget), 'UniformOutput', 0))';
                    lfp2_tmp = cell2mat(cellfun(@transpose, lfp2(indstoget), 'UniformOutput', 0))';
                    
                    [C, t, f] = lt_neural_Cohere_Recalcsub(thingstodo, lfp1_tmp, lfp2_tmp, cohversion, ...
                        ntapers, movingwin, tw, t_LFP);
                    
                    if isempty(t)
                        t = t_LFP;
                    end
                    
                    assert(isempty(CohMean_Base{indthis_out}));
                    CohMean_Base{indthis_out} = C;
                    
                    if any(strcmp(thingstodo, 'lfpxcorr'))
                        [xcorr_max, xcorr_lag, freqsthis] = ...
                            lt_neural_LFP_RecalcCoh_sub0(filtdatstruct, ...
                            chan1, chan2, indstoget, ...
                            lfpxcorr_twindtoget, lfpxcorr_ffwindtoget, onlyOnePeriod, fs, ...
                            usecorrcoeff);
                        
                        tmp = nan(size(xcorr_max,2), size(xcorr_max,1), 2); % trials x fvals x 2[c and lag]
                        tmp(:,:, 1) = xcorr_max';
                        tmp(:,:, 2) = xcorr_lag';
                        LFPXCORR_Base_tr_f_cclag{indthis_out} = tmp;
                        LFPXCORR_freqsall{indthis_out} = freqsthis;
                    end
                    
                    
                    
                    
                    % ##########################  WN
                    indstoget = inds_WN;
                    
                    lfp1_tmp = cell2mat(cellfun(@transpose, lfp1(indstoget), 'UniformOutput', 0))';
                    lfp2_tmp = cell2mat(cellfun(@transpose, lfp2(indstoget), 'UniformOutput', 0))';
                    
                    [C, t, f] = lt_neural_Cohere_Recalcsub(thingstodo, lfp1_tmp, lfp2_tmp, cohversion, ...
                        ntapers, movingwin, tw, t_LFP);
                    if isempty(t)
                        t = t_LFP;
                    end
                    
                    assert(isempty(CohMean_WN{indthis_out}));
                    CohMean_WN{indthis_out} = C;
                    
                    if any(strcmp(thingstodo, 'lfpxcorr'))
                        [xcorr_max, xcorr_lag, freqsthis] = ...
                            lt_neural_LFP_RecalcCoh_sub0(filtdatstruct, ...
                            chan1, chan2, indstoget, ...
                            lfpxcorr_twindtoget, lfpxcorr_ffwindtoget, onlyOnePeriod, fs, ...
                            usecorrcoeff);
                        
                        tmp = nan(size(xcorr_max,2), size(xcorr_max,1), 2); % trials x fvals x 2[c and lag]
                        tmp(:,:, 1) = xcorr_max';
                        tmp(:,:, 2) = xcorr_lag';
                        LFPXCORR_WN_tr_f_cclag{indthis_out} = tmp;
                    end

                    
                    
                    
                    % -=-- sanity check, compare previuosly extracted to
                    % current.
                    if (0)
                        lt_figure; hold on;
                        lt_subplot(4,2,1); hold on;
                        title('old [base]')
                        lt_neural_Coher_Plot(OUTSTRUCT.CohMean_Base{indthis_out}, ...
                            PARAMS.tbins, PARAMS.ffbins, 1, '', [0.2 0.8], 0);
                        
                        lt_subplot(4,2,2); hold on;
                        title('new [base]')
                        lt_neural_Coher_Plot(CohMean_Base{indthis_out}, ...
                            PARAMS.tbins, PARAMS.ffbins, 1, '', [0.2 0.8], 0);
                        
                        lt_subplot(4,2,3); hold on;
                        title('old [wn]')
                        lt_neural_Coher_Plot(OUTSTRUCT.CohMean_WN{indthis_out}, ...
                            PARAMS.tbins, PARAMS.ffbins, 1, '', [0.2 0.8], 0);
                        
                        lt_subplot(4,2,4); hold on;
                        title('new [wn]')
                        lt_neural_Coher_Plot(CohMean_WN{indthis_out}, ...
                            PARAMS.tbins, PARAMS.ffbins, 1, '', [0.2 0.8], 0);
                        
                        lt_subplot(4,2,5); hold on;
                        title('base');
                        plot(OUTSTRUCT.CohMean_Base{indthis_out}(:), ...
                            CohMean_Base{indthis_out}(:), 'ok');
                        
                        lt_subplot(4,2,6); hold on;
                        title('wn');
                        plot(OUTSTRUCT.CohMean_WN{indthis_out}(:), ...
                            CohMean_WN{indthis_out}(:), 'ok');
                    end
                    % --- collect t, f
                    Tall = [Tall; t];
                    Fall = [Fall; f];
                    
                end
                
            end
        end
    end
end

if all(~cellfun(@isempty, CohMean_Base))
    OUTSTRUCT.CohMean_Base = CohMean_Base;
    OUTSTRUCT.CohMean_WN = CohMean_WN;
% UPDATE PARAMS
PARAMS.tbins = Tall(1,:);
PARAMS.ffbins = Fall(1,:);
end


% === plot sampel size
lt_figure; hold on;
xlabel('base');
ylabel('wn');
plot(SampleSizes(:,1), SampleSizes(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b', []);
