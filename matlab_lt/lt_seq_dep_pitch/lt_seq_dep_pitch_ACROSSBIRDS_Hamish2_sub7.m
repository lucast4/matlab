function DATSTRUCT = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub7(DATSTRUCT, Indstoplot, ...
    NrendsORIG, usemedian, normwiggle, doBase)

% doBase=1; % default = 1; if 0 then does training epioch
if ~exist('doBase', 'var')
    doBase = 1;
end

%%
% Indstoplot = find(All_Issame==1 | All_Issame==0); % uses all syllables
% NrendsORIG = 25; % renditions to take from edges (will minimize if not enough trials
% usemedian = 1; % if 0, then uses means. if 1, then uses MAD and medians
% normwiggle = 1; % if 1, then normalized wiggle to the wiggle for MUSC (matched sample size)

%%

%% extract variables

fnames = fieldnames(DATSTRUCT)';
for fn=fnames
    fn = fn{1};
    eval([fn ' = DATSTRUCT.(fn);'])
end

%% baseline or trainnig?
if doBase==1
    fnamethis = 'All_PitchCont_BASE_PBS';
    fnamethis_MUSC = 'All_PitchCont_BASE_MUSC';
elseif doBase==0
    % - then do training
    fnamethis = 'All_PitchCont_WN_PBS';
    fnamethis_MUSC = 'All_PitchCont_WN_MUSC';
end

%% =========== COLLECT THINGS
WiggleAll = cell(length(Indstoplot), 2); % syls x [lo, hi]
WiggleVsMeanFF = nan(length(Indstoplot), 1);
BaseBiasAll = nan(length(Indstoplot),1);
TwindDurAll = nan(length(Indstoplot),1);
WiggleOverMusc = nan(length(Indstoplot), 1);
Wiggle = nan(length(Indstoplot), 1);


for i=1:length(Indstoplot)
    
    indthis = Indstoplot(i);
    
    %  ########################## COLLECT DAY SUMMARIES AND THEN AVERAGE
    numdays = length(DATSTRUCT.(fnamethis)(indthis).All_PCmat);
    
    rho_allday = []; % initiate, will collect across days then take average.
    wiggle_all_allday = [];
    wiggle_hi_allday = [];
    wiggle_lo_allday = [];
    afpbias_allday = [];
    datagood = 1;
    
    for dd = 1:numdays
        
        % ============================ BASELINE, PBS
        pcmat = DATSTRUCT.(fnamethis)(indthis).All_PCmat{dd};
        ffmat = DATSTRUCT.(fnamethis)(indthis).All_ffvals{dd};
        tthis = DATSTRUCT.(fnamethis)(indthis).All_tbins{dd};
        twind = DATSTRUCT.(fnamethis)(indthis).All_twind(dd,:);
        
        pcmat_MUSC = DATSTRUCT.(fnamethis_MUSC)(indthis).All_PCmat{dd};
        
        if any(isnan(twind))
            % then no time window defined for this syl, shoudl skuip
            datagood = 0;
            break
        end
        
        %         twind(1) = twind(1)+40;
        %         twind(2) = twind(2)-40
        
        % ============================= COLLECT PC AT EXTREMES
        % ------------- minimum sample size
        if length(ffmat)<1*NrendsORIG
            continue
        end
        if length(ffmat)<2*NrendsORIG
            % then make nrends smaller
            Nrends = floor(length(ffmat)/2);
        else
            Nrends = NrendsORIG;
        end
        
        % ------- match to musc, if desired
        if normwiggle ==1
            
            % --- cut down sample size to match that of musc (if musc few)
            n_musc = size(pcmat_MUSC,1);
            Nrends = min([Nrends, n_musc]);
            pcmat_MUSC = pcmat_MUSC(randperm(size(pcmat_MUSC,1), Nrends), :);
            
            % -- normalize
            pcmat_MUSC = pcmat_MUSC - mean(pcmat_MUSC, 1); % subtract x-trial mean
            pcmat_MUSC = pcmat_MUSC - mean(pcmat_MUSC(:, twind(1):twind(end)), 2); % subtract within-trial means
            
            % -- collect wiggle for MUSC trials, then take avearge across
            % trials.
            if usemedian==1
                wiggle_musc =  mad(pcmat_MUSC(:, twind(1):twind(end)), [], 2);
                wiggle_musc = median(wiggle_musc);
            else
                wiggle_musc = std(pcmat_MUSC(:, twind(1):twind(end)), [], 2);
                wiggle_musc = mean(wiggle_musc);
            end
        end
        
        
        
        % ========== 1) NORMALIZE ALL PCS: 1) SUBTRACT X-TRIAL MEAN, 2) SUBTRACT WITHIN TRIAL MEAN
        % ---- subtract x-trial mean
        pcmat_norm = pcmat - mean(pcmat, 1);
        
        % ---- subtract within trial means
        pcmat_norm = pcmat_norm - mean(pcmat_norm(:, twind(1):twind(end)), 2); % subtract within trial mean
        
        % ---- COLLECT WIGGLES FOR EACH TRIAL
        if usemedian==1
            wiggles_all = mad(pcmat_norm(:, twind(1):twind(end)), [], 2);
        else
            wiggles_all = std(pcmat_norm(:, twind(1):twind(end)), [], 2);
        end
        
        % ============== COLLECT ALL WIGGLES
        if usemedian==1
            wiggle = median(wiggles_all);
        else
            wiggle = mean(wiggles_all);
        end
        if normwiggle==1
            wiggle = wiggle./wiggle_musc;
        end
        wiggle_all_allday = [wiggle_all_allday wiggle];
        
        
        % =============== 2) COLLECT HI AND LO DATASETS
        [~, indsort] = sort(ffmat);
        
        % ---- hi pitch
        if usemedian==1
            wiggle = median(wiggles_all(indsort(end-Nrends+1:end), :));
        else
            wiggle = mean(wiggles_all(indsort(end-Nrends+1:end), :));
        end
        if normwiggle==1
            wiggle = wiggle./wiggle_musc;
        end
        wiggle_hi_allday = [wiggle_hi_allday wiggle];
        %         WiggleAll{i,2} = wiggle;
        
        % ---- lo pitch
        if usemedian ==1
            wiggle = median(wiggles_all(indsort(1:Nrends), :));
        else
            wiggle = mean(wiggles_all(indsort(1:Nrends), :));
        end
        if normwiggle==1
            wiggle = wiggle./wiggle_musc;
        end
        wiggle_lo_allday = [wiggle_lo_allday wiggle];
        %         WiggleAll{i,1} = wiggle;
        
        
        % ================= 3) correlation coeff between ff and wiggle
        if usemedian==1
            rho = corr(tiedrank(wiggles_all), tiedrank(ffmat)); % rank corr
        else
            rho = corr(wiggles_all, ffmat);
        end
        rho_allday = [rho_allday rho];
        %         WiggleVsMeanFF(i) = rho;
        
        
        
        % ######################### BASELINE AFP BIAS?
        pcmat_PBS = All_PitchCont_BASE_PBS(indthis).All_PCmat{dd};
        pcmat_MUSC = All_PitchCont_BASE_MUSC(indthis).All_PCmat{dd};
        
        ffPBS = mean(pcmat_PBS(:, twind(1):twind(end)),2);
        ffMUSC = mean(pcmat_MUSC(:, twind(1):twind(end)),2);
        
        afpbias = mean(ffPBS) - mean(ffMUSC);
        %         afpbias = mean(ffPBS) - mean(ffMUSC);
        
        % --- collect across days
        afpbias_allday = [afpbias_allday afpbias];
        %         BaseBiasAll(i) = afpbias;
        
    end
    
    % ========== take average of day stats
    if datagood==1
        WiggleVsMeanFF(i) = mean(rho_allday);
        BaseBiasAll(i) = mean(afpbias_allday);
        WiggleAll{i,1} = mean(wiggle_lo_allday);
        WiggleAll{i,2} = mean(wiggle_hi_allday);
        Wiggle(i) = mean(wiggle_all_allday);
        TwindDurAll(i) = tthis(twind(2)) - tthis(twind(1));
    end
    
    if isempty(WiggleAll{i,2})
        WiggleOverMusc(i) = nan;
    else
        WiggleOverMusc(i) = WiggleAll{i,2}>1;
    end
end



% ################################################################ OUTPUTS
DATSTRUCT.BaseBiasAll = BaseBiasAll;
DATSTRUCT.WiggleVsMeanFF = WiggleVsMeanFF;
DATSTRUCT.WiggleAll = WiggleAll;
DATSTRUCT.Wiggle = Wiggle;
DATSTRUCT.TwindDurAll = TwindDurAll;
DATSTRUCT.WiggleOverMusc = WiggleOverMusc;

%

disp('DONE ! ...');
