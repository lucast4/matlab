%% ================ extract ngrams and fr
close all;

Params.LearnKeepOnlyBase = 1;
Params.strtype = 'xaaa';
Params.Nmin = 6; % num trials minimum
Params.alignsyl = 3;

Params.regexpr.motifpredur = 0.2;
Params.regexpr.motifpostdur = 0.1;
Params.regexpr.alignOnset = 1;
Params.regexpr.preAndPostDurRelSameTimept=1;
Params.regexpr.RemoveIfTooLongGapDur = 1;

saveON =1;

NGRAMSTRUCT = lt_neural_NGRAMS_Extract(SummaryStruct, Params, saveON);

%% #################################### EXTRACTION/PREPROCESSING
%% ================= EXTRACT NGRAMSTRUCT FROM SAVED DATA 
close all;
dirname = 'xaa_30Apr2018_2342';
window_prem = [-0.025 0.025]; % relative to syl onset
Nshuffs = 1; % for negative control (corr analysis)
doSqrtTransform = 1;
use_dPrime = 0; % if 1, then gets mean dPrime instead of mean abs FR diff
nshufftmp = 2; % for zscoring fr diff.
DoDecode = 1; % to get decoder performances; IMPORTNAT: IF 1, then overwrites 
% FRdiff values.

[OUTSTRUCT, SummaryStruct, Params] = lt_neural_NGRAMS_Compile(dirname, ...
    window_prem, Nshuffs, doSqrtTransform, use_dPrime, nshufftmp, DoDecode);
Params.window_prem = window_prem;
Params.doSqrtTransform = doSqrtTransform;
Params.use_dPrime = use_dPrime;

% ========================== FOR COMPATIBILITY WITH OLD CODE, EXTRACT ALL
% FIELDS
fnamesthis = fieldnames(OUTSTRUCT);
for j=1:length(fnamesthis)
   eval([fnamesthis{j} ' = OUTSTRUCT.' fnamesthis{j} ';']); 
end

%% ================= get bregion for each datapt
All_Bregion = cell(length(OUTSTRUCT.All_birdnum),1);
for i=1:length(OUTSTRUCT.All_birdnum)
   bnum = OUTSTRUCT.All_birdnum(i);
   neur = OUTSTRUCT.All_neurnum(i);
   disp(num2str(i));
   All_Bregion{i} = SummaryStruct.birds(bnum).neurons(neur).NOTE_Location;
end
OUTSTRUCT.All_Bregion = All_Bregion;

%% ================= get all pairwise distances during premotor window
if (0) % OLD VERSION --- this works with NGRAMSTRUCT. new version does not since 
    % filesize too large.
    
    % if ever want to use this version, need to run the script in here.
    % this extracts summary arrays.
lt_neural_NGRAMS_GetDatOldVersion;

end
%% ================== convert motifpair types to groups

if strcmp(Params.strtype, 'xaa')
    PairTypesInOrder = {...,
        '0  0  1', ...
        '1  0  0', ...
        '1  0  1', ...
        '0  1  0', ...
        '0  1  1', ...
        '1  1  0', ...
        '1  1  1'}; % in order to be plotted
elseif strcmp(Params.strtype, 'xaaa')
    PairTypesInOrder = {...,
        '0  0  0  1', ...
        '1  0  0  0', ...
        '1  0  0  1', ...
        '0  1  0  0', ...
        '0  1  0  1', ...
        '1  1  0  0', ...
        '1  1  0  1', ...
        '0  0  1  1', ...
        '0  0  1  0', ...
        '0  1  1  0', ...
        '0  1  1  1', ...
        '1  1  1  0', ...
        '1  0  1  0', ...
        '1  0  1  1', ...
        '1  1  1  1', ...
        }; % in order to be plotted
end
All_diffsyl_string = num2str(double(All_diffsyl_logical));
All_diffsyl_string = mat2cell(All_diffsyl_string, ones(size(All_diffsyl_logical,1),1));

[~, All_diffsyl_PairType] = ismember(All_diffsyl_string, PairTypesInOrder);
assert(all(strcmp(All_diffsyl_string, PairTypesInOrder(All_diffsyl_PairType)')), 'asdfas');

% ================= PUT INTO STRUCT
OUTSTRUCT.PairTypesInOrder = PairTypesInOrder;
OUTSTRUCT.All_diffsyl_PairType = All_diffsyl_PairType;


%% ============ [FIGURE OUT MISLABELED SYLLABLES]

ncases = length(OUTSTRUCT.All_birdnum);
OUTSTRUCT.All_BadSyls = nan(ncases,1);
for i=1:ncases
    disp(num2str(i));
    % -- birdname
    birdname = SummaryStruct.birds(OUTSTRUCT.All_birdnum(i)).birdname;    
    badstrings = lt_neural_NGRAMS_BadSyls(birdname, Params);
    
    % -- if either of then grams contains any of the bad strings, then
    % mark as bad
    ngramstrings = OUTSTRUCT.All_ngramstring_inorder(i,:);
    
%     if strcmp(birdname, 'pu26y2')
%         keyboard
%     end
    anybad = [];
    for j=1:2
        tmp = regexp(ngramstrings{j}, badstrings);
        anybad = any([anybad ~cellfun(@isempty, tmp)]);
    end
    
    OUTSTRUCT.All_BadSyls(i) = anybad;
end


% =========== for each bird, list all motif pairs that are bad
if (0)
    numbirds = max(OUTSTRUCT.All_birdnum);
for j=1:numbirds
   inds = OUTSTRUCT.All_birdnum==j & OUTSTRUCT.All_BadSyls==1;
   motifpairs = OUTSTRUCT.All_ngramstring_inorder(inds,:);
  birdname = SummaryStruct.birds(j).birdname;
  
   disp([' ==================== ' birdname]);
   
   motifpairstrings = cell(size(motifpairs,1), 1);
   for jj=1:size(motifpairs,1)
       motifpairstrings{jj} = [motifpairs{jj,1} '-' motifpairs{jj,2}];
   end
   disp(unique(motifpairstrings));
   disp(' ---- ignoring any with "x"');
   tmp = unique(motifpairstrings);
   indsnox = cellfun(@isempty, regexp(tmp, 'x'));
   disp(tmp(indsnox));
end

end

%% ============ REMOVE BAD LABEL PAIRS FROM DATASET

% ==== save original outsturct
OUTSTRUCT_orig = OUTSTRUCT;

% ==== save this for later
PairTypesInOrder = OUTSTRUCT.PairTypesInOrder;

% ==== remove this field, since is smaller vector.
OUTSTRUCT = rmfield(OUTSTRUCT, 'PairTypesInOrder');

% ==== go thru all fields and only keep the good syls
indstokeep = ~OUTSTRUCT.All_BadSyls;
fnames = fieldnames(OUTSTRUCT);

for j=1:length(fnames)
    
    ytmp = OUTSTRUCT.(fnames{j});
    
    try
        ytmp = ytmp(indstokeep, :);
        OUTSTRUCT.(fnames{j}) = ytmp;
    catch err
        try ytmp = ytmp(:, indstokeep);
            OUTSTRUCT.(fnames{j}) = ytmp';
        catch err
            disp('why error?');
        end
    end
end

% === put pairtypes back in
OUTSTRUCT.PairTypesInOrder = PairTypesInOrder;

%% =========== [RE-EXTRACT, EQUALIZING SAMPLE SIZE]
close all;
measure_to_recalc = 'absfrdiff';
if strcmp(Params.strtype, 'xaa')
PairTypesToCompare = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
 PairTypesToCompare = {...
    '1  1  1  1', ... % xaxis
    '1  0  0  0'}; % yaxis
end   
nshufftmp = 2;
DoDecode =1; % IMPORTANT: if this is 1, then uses decode and overwrites FR diff stuff
OUTSTRUCT = lt_neural_NGRAMS_ReSample(OUTSTRUCT, SummaryStruct, Params, ...
    measure_to_recalc, PairTypesToCompare, nshufftmp, DoDecode);

%% ============= save outstruct
savesuffix = 'DecodeWithSqrtTransform';

fname = ['/bluejay5/lucas/analyses/neural/NGRAMS/' Params.dirname '/OUTSTRUCT_' savesuffix '.mat'];
save(fname, 'OUTSTRUCT');

fname = ['/bluejay5/lucas/analyses/neural/NGRAMS/' Params.dirname '/Params_' savesuffix '.mat'];
save(fname, 'Params');

%% ============ [DOWNSAMPLING]
% QUESTION: RA vs. LMAN, what if downsample RA so that equalize positive
% control values between RA and LMAN?

% === for each neuron, regardless of brain region, get values as function
% of sample size.

measure_to_recalc = 'absfrdiff';
downfactorlist = [0.025 0.05 0.1 0.2 0.5 0.75 1]; % fractin of samples between N and Nmin.
[OUTSTRUCT_subsamp, downfactorlist] = lt_neural_NGRAMS_Downsample(OUTSTRUCT, SummaryStruct, Params, ...
    measure_to_recalc, downfactorlist);
Params.downfactorlist = downfactorlist;

%% ============= [DOWNSAMPLE] === plot each neuron, fr diff 
% as function of downsample.

close all;

PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
Indpaircomp = find(ismember(OUTSTRUCT.PairTypesInOrder, PairtypesToplot));

% ============
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);


for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    for ii=1:maxneur
        
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        if ~any(inds)
            continue
        end
        
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        
        % =====================
        PairTypes = OUTSTRUCT.All_diffsyl_PairType(inds);
        Ydat = OUTSTRUCT_subsamp.FRdiffDAT(inds,:);
        Yshuff = OUTSTRUCT_subsamp.FRdiffShuff(inds,:);
        Nmean = OUTSTRUCT_subsamp.Nboth(inds, :, :);
        Nmean = squeeze(mean(Nmean, 2));
        
        % --------------------------- first pairtype (x)
        indtype = 1;
        
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = Ydat(indtmp,:);
        yshuff = Yshuff(indtmp, :);
%         x = Nmean(indtmp,:); % sample size, mean between 2 motifs
        x = Params.downfactorlist;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-n' num2str(ii) '-' bregion]);
%         xlabel('subsample N (mean 2 motifs)');
        xlabel('downfactor');
        ylabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)}]);
        % -- dat
        if sum(indtmp)<100
        plot(x, ydat, '-', 'Color', [0.7 0.7 0.7]);
        plot(x, yshuff, '-', 'Color', [0.7 0.2 0.2]);
        end
        % -- means
        plot(x, mean(yshuff,1), '-r', 'LineWidth', 2);
        plot(x, mean(ydat,1), '-k', 'LineWidth', 2);
        
        
        % --------------------------- second pairtype (x)
        indtype = 2;
        
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = Ydat(indtmp,:);
        yshuff = Yshuff(indtmp, :);
%         x = Nmean(indtmp,:); % sample size, mean between 2 motifs
        x = Params.downfactorlist;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-n' num2str(ii)]);
%         xlabel('subsample N (mean 2 motifs)');
        xlabel('downfactor');
        ylabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)}]);
        
        % -- dat
        if sum(indtmp)<100
        plot(x, ydat, '-', 'Color', [0.7 0.7 0.7]);
        plot(x, yshuff, '-', 'Color', [0.7 0.2 0.2]);
        end
        % -- means
        plot(x, mean(yshuff,1), '-r', 'LineWidth', 2);
        plot(x, mean(ydat,1), '-k', 'LineWidth', 2);
        
        
        % ------------------------ COMBINE (EACH TYPE, SUBTRACT GLOBAL NEG)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-n' num2str(ii)]);
%         xlabel('subsample N (mean 2 motifs)');
        xlabel('downfactor');
        ylabel(['frdiff']);
        
        Yall = {}; % dat, dat, shuff
        
        % -- type 1 (dat)
        indtype = 1;
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = mean(Ydat(indtmp,:),1);
        plot(x, ydat, '--r');
        lt_plot_text(x(1), ydat(1), OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)});
                
        Yall{indtype} = ydat;
        
        % -- type 1 (dat)
        indtype = 2;
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = mean(Ydat(indtmp,:),1);
        plot(x, ydat, '--b');
        lt_plot_text(x(1), ydat(1), OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)});
                
        Yall{indtype} = ydat;
        
        % -- negative
        indtmp = ismember(PairTypes, Indpaircomp);
        yshuff = mean(Yshuff(indtmp,:),1);
        
        plot(x, yshuff, '--k')
        
        % ============= OVERLAY DIFFERENCES
        plot(x, Yall{1}-yshuff, '-r');
        plot(x, Yall{2}-yshuff, '-b');
        lt_plot_zeroline;
    end
end


%% =========== [DOWNSAMPLE] - replace by certain factor of downsample
% then redo scatter plot analysis
dsampvalue_take = 0.05;
dsampbregion = {'RA'}; % if empty, downsamples none. e.g. {'LMAN', 'RA'};

OUTSTRUCT_tmp = OUTSTRUCT;

dsampind = find(Params.downfactorlist ==  dsampvalue_take);
assert(length(dsampind)==1, 'this dsamp value doesnt exist...');

% ================= REPLACE DATA WITH DOWNSAMPLED DATA
indstodo = ismember(OUTSTRUCT_tmp.All_Bregion, dsampbregion);

OUTSTRUCT_tmp.All_AbsFRdiff(indstodo) = ...
    OUTSTRUCT_subsamp.FRdiffDAT(indstodo, dsampind);
OUTSTRUCT_tmp.All_AbsFRdiff_NEG(indstodo) = ...
    OUTSTRUCT_subsamp.FRdiffShuff(indstodo, dsampind);
OUTSTRUCT_tmp.All_AbsFRdiff_Zrelshuff(indstodo) = ...
    OUTSTRUCT_subsamp.FRdiff_Z(indstodo, dsampind);
OUTSTRUCT_tmp.All_N(indstodo, :) = ...
    OUTSTRUCT_subsamp.Nboth(indstodo, :, dsampind);

% ================= RUN SCATTERPLOT
close all;

plottype = 'absfrdiff_typediff'; % oneminusrho or absfrdiff or absfrdiff_globZ or absfrdiff_typediff
usemedian = 0; % only workds for absfrdiff, absfrdiff_globZ or absfrdiff_typediff
plotON=0; % only works for absfrdiff
if strcmp(Params.strtype, 'xaa')
PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
 PairtypesToplot = {...
    '1  1  1  1', ... % xaxis
    '1  0  0  0'}; % yaxis
end   
plotRawGood = 0; % histogram for pos, negative, dat (only works for plottype = absfrdiff_globZ);

% ----------- params for one minus rho, specifically
dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
sanitycheckuseneg = 0; % uses negative control data instead of data
removeBadSyls = 1; % i.e. badly labeled...

[AllPairs_Means, AllPairs_Birdnum, AllPairs_Bregions] = ...
    lt_neural_NGRAMS_PlotScatter(OUTSTRUCT_tmp, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg, plotRawGood, usemedian, ...
    removeBadSyls);




%% ####################################################
%% ============ [DIAGNOSTICS]
close all;
dispNgramStrings = 0; % then for 10% of neurons will disp.
plotRawSampSize = 0; % then plots number of trials.
compareNtoFRdist = 1; % MANY PLOTS - compares samp size to fr dist. neuron by neuron.
lt_neural_NGRAMS_DIAGNOSTIC(OUTSTRUCT, SummaryStruct, Params, dispNgramStrings, ...
    plotRawSampSize, compareNtoFRdist);


%% ============ [PLOT] separate by pair type, and average across units
close all;

plottype = 'absfrdiff'; % oneminusrho or absfrdiff
plotON=1; % raw plots? only works for absfrdiff
dosubtractcontrol = 0; % then subtracts negative control before plotting

lt_neural_NGRAMS_PlotByPairtype(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    dosubtractcontrol)


%% =========== [SCATTER PLOT] for each bird, plot scatter of two pairtypes

% ALSO PLOTS HISTOGRAM
close all;

plottype = 'absfrdiff_globZ'; % oneminusrho or absfrdiff or absfrdiff_globZ or absfrdiff_typediff
usemedian = 0; % only workds for absfrdiff, absfrdiff_globZ or absfrdiff_typediff
plotON=0; % only works for absfrdiff
if strcmp(Params.strtype, 'xaa')
PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '0  1  1'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
 PairtypesToplot = {...
    '1  1  1  1', ... % xaxis
    '1  0  0  0'}; % yaxis
end   
Nmin = 2; % number of pairs in class. [DOESNT WORK YET]
plotRawGood = 1; % histogram for pos, negative, dat (only works for plottype = absfrdiff_globZ);

% ----------- params for one minus rho, specifically
dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
sanitycheckuseneg = 0; % uses negative control data instead of data
removeBadSyls = 1; % i.e. badly labeled...

[AllPairs_Means, AllPairs_Birdnum, AllPairs_Bregions] = ...
    lt_neural_NGRAMS_PlotScatter(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg, plotRawGood, usemedian, ...
    removeBadSyls);


%%  ######################## REGRESSION MODELING


lt_neural_NGRAMS_Regression;


%% ########################### PLOT EXAMPLE FR TRACES
close all;
birdtoplot = 'pu69wh78';
neurtoplot = 2;
ngramstoplot = {}; % leave empty to plot random one
pairtypetoplot = '1  1  1';

 lt_neural_NGRAMS_PlotEgPair(OUTSTRUCT, SummaryStruct, Params, ...
    birdtoplot, neurtoplot, ngramstoplot, pairtypetoplot);

%% =========== [DIAGNOSTIC] - plot motif pair strings...
close all;

plottype = 'absfrdiff_globZ'; % oneminusrho or absfrdiff or absfrdiff_globZ or absfrdiff_typediff
if strcmp(Params.strtype, 'xaa')
PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
 PairtypesToplot = {...
    '1  1  1  1', ... % xaxis
    '1  0  0  0'}; % yaxis
end   

% -------------- WHICH TO PLOT?
birdtoplot = 'pu69wh78';
neurtoplot = 14;

lt_neural_NGRAMS_PlotMotStr(OUTSTRUCT, SummaryStruct, plottype, ...
    PairtypesToplot, birdtoplot, neurtoplot);

%% ============ [HISTOGRAMS] FOR A GIVEN NEURON, PLOT HISTOGRAMS
close all;
PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
assert(strcmp(plottype, 'absfrdiff_globZ'), 'asfdsa')
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);

for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    for ii=1:maxneur
        
        inds = All_birdnum==i & All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;


%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%             title([birdname '-n' num2str(ii) '[' bregion ']']);

        
        for k =1:numpairtypes
            
            pairtypethis = PairtypesToplot{k};
            pairtypethis = find(strcmp(PairTypesInOrder, pairtypethis));
            
            inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
            
            
            % ########################################## COLLECT DATA
            y = [];
            if strcmp(plottype, 'absfrdiff')
                % then mean abs diff in FR, normalized within each
                % pairtype
                y = All_AbsFRdiff_Zrelshuff(inds);
                
            elseif strcmp(plottype, 'absfrdiff_globZ')
                % =========== version with both dat and pos compared to
                % distribution of negative controls
                
                % ----- V1 - NEG = combined from the pairs being analyzed
                if useGlobalNeg==0 % old version, limited to just the desired pairs. but I think
                    % is better to use all pairs ...
                    [~, pairtypes_neg] = intersect(PairTypesInOrder, PairtypesToplot);
                    indsneg = All_birdnum==i & All_neurnum==ii & ...
                        ismember(All_diffsyl_PairType, pairtypes_neg);
                elseif useGlobalNeg==1
                    indsneg = All_birdnum==i & All_neurnum==ii;
                end
                
                negmean = mean(All_AbsFRdiff_NEG(indsneg));
                negstd = std(All_AbsFRdiff_NEG(indsneg));
                
                y = (All_AbsFRdiff(inds) - negmean)./negstd;
                
            elseif strcmp(plottype, 'absfrdiff_typediff')
                
                % ----- V2 - each type compaired to mean of its own neg
                % (not zscored)
                y = All_AbsFRdiff(inds) - mean(All_AbsFRdiff_NEG(inds));
                
            end
            
            
            % ###################################### PLOT HISTOGRAMS?
            if plotON==1
                % =============== plot histograms
                lt_plot_histogram(y, '', 1, 1, '', 1, pcolors{k});
            end
            
            
            % #################### collect mean to then do scatterplot
            ymean = mean(y);
            allmeans(k) = ymean;
        end
    end
end


%% =========== [CHECK] Compare resampled to old data, separated by pairtype
% === one plot per neuron, grouping by pairtype, comparing oold vs. new
% === only plots fr diff (for both data and shuffle).
close all; 
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);

figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    for ii=1:maxneur
   
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        if ~any(inds)
            continue
        end
        
        % ===================== 
        PairTypes = OUTSTRUCT.All_diffsyl_PairType(inds);
        
        Y_old = OUTSTRUCT.All_AbsFRdiff_ORIG(inds);
        Yneg_old = OUTSTRUCT.All_AbsFRdiff_NEG_ORIG(inds);
        
        Y_new = OUTSTRUCT.All_AbsFRdiff(inds);
        Yneg_new = OUTSTRUCT.All_AbsFRdiff_NEG(inds);
        
        % ==================================== PLOTS [DAT]
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-' num2str(ii)]);
        xlabel('old(gr) - new(bk)');
        ylabel('fr dist');
        
        % ==== plot old
        plot(PairTypes-0.2, Y_old, 'x', 'Color', [0.7 0.7 0.7]);
        [ymean, ystd] = grpstats(Y_old, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x-0.1, ymean, {'Errors', ystd, 'Color', [0.7 0.7 0.7]});
        
        % ==== plot new
        plot(PairTypes+0.2, Y_new, 'xk');
        [ymean, ystd] = grpstats(Y_new, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x+0.3, ymean, {'Errors', ystd, 'Color', 'k'});

        % ----
        set(gca, 'XTick', 1:max(PairTypes));
        
        % ==================================== PLOTS [NEG SHUFF]
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['[SHUFF]' birdname ',' num2str(i) '-' num2str(ii)]);
        xlabel('old(gr) - new(bk)');
        ylabel('fr dist');
        
        % ==== plot old
        plot(PairTypes-0.2, Yneg_old, 'x', 'Color', [0.7 0.7 0.7]);
        [ymean, ystd] = grpstats(Yneg_old, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x-0.1, ymean, {'Errors', ystd, 'Color', [0.7 0.7 0.7]});
        
        % ==== plot new
        plot(PairTypes+0.2, Yneg_new, 'xk');
        [ymean, ystd] = grpstats(Yneg_new, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x+0.3, ymean, {'Errors', ystd, 'Color', 'k'});

        % ----
        axis tight;
        set(gca, 'XTick', 1:max(PairTypes));
        lt_plot_zeroline;
    end
    
    if mod(i,5)==0
        disp('PRESS ANYTHING TO GO TO NEXT SET OF 5 BIRDS')
        pause
        close all;
    end
end
    
    
%% ==== [CHECK] compare sample size and negative control decoding
% ===== for two pairtypes, does scatterplot, with each datapoint one
% neuron.
% -== one plot per animal.
% === can compare things like sample size.

PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
Indpaircomp = find(ismember(OUTSTRUCT.PairTypesInOrder, PairtypesToplot));

plottype = 'N_new';
% ---- options:
% frdiff_shuff_new
% frdiff_shuff_old
% N_old
% N_new

% ============
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);

figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    % ================== one subplot per bird
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' plottype ']' birdname ',' num2str(i)]);
    xlabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(1)}]);
    ylabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(2)}]);
    
    YYall = [];
    YYstdall = [];
    for ii=1:maxneur
        
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        if ~any(inds)
            continue
        end
        
        % =====================
        PairTypes = OUTSTRUCT.All_diffsyl_PairType(inds);
        
        if strcmp(plottype, 'frdiff_shuff_new')
            Ythis = OUTSTRUCT.All_AbsFRdiff_NEG(inds);
        elseif strcmp(plottype, 'frdiff_shuff_old')
            Ythis = OUTSTRUCT.All_AbsFRdiff_NEG_ORIG(inds);
        elseif strcmp(plottype, 'N_old')
            Ythis = OUTSTRUCT.All_N_ORIG(inds, :);
            Ythis = mean(Ythis,2); % combine across ngrams
        elseif strcmp(plottype, 'N_new')
            Ythis = OUTSTRUCT.All_N(inds, :);
            Ythis = mean(Ythis,2); % combine across ngrams
        end
%         
%         Y_old = OUTSTRUCT.All_AbsFRdiff_ORIG(inds);
%         Yneg_old = OUTSTRUCT.All_AbsFRdiff_NEG_ORIG(inds);
%         
%         Y_new = OUTSTRUCT.All_AbsFRdiff(inds);
%         Yneg_new = OUTSTRUCT.All_AbsFRdiff_NEG(inds);
%         
%         N_old = OUTSTRUCT.All_N(inds,:);
%         N_new = OUTSTRUCT.All_N_ORIG(inds,:);
%         
        
        % ========================= PLOTS [neg shuff distribition];
        
        % ========== NEW
        Y_all = {};
        Ystd_all = {};
        
        % --------------------------- first pairtype (x)
        indtmp = PairTypes == Indpaircomp(1);
        y = Ythis(indtmp);
        
        % -- collect
        Y_all{1} = mean(y);
        Ystd_all{1} = std(y);
        
        % --------------------------- second pairtype (x)
        indtmp = PairTypes == Indpaircomp(2);
        y = Ythis(indtmp);
        
        % -- collect
        Y_all{2} = mean(y);
        Ystd_all{2} = std(y);
        
        % ======================
        YYall = [YYall; [Y_all{1} Y_all{2}]];
        YYstdall = [YYstdall ; [Ystd_all{1} Ystd_all{2}]];

        % ==================== PLOT
%         lt_plot(Y_all{1}, Y_all{2});
        line([Y_all{1} Y_all{1}], [Y_all{2}-Ystd_all{2} Y_all{2}+Ystd_all{2}], 'Color', 'k');
        line([Y_all{1}-Ystd_all{1} Y_all{1}+Ystd_all{1}], [Y_all{2} Y_all{2}], 'Color', 'k');
    end
    
    % =========== PLOT
%     lt_plot_45degScatter(YYall(:,1), YYall(:,2), 'k', 1);
    lt_plot(YYall(:,1), YYall(:,2), {'Color', 'k'});
    YLIM = ylim;
    XLIM = xlim;
    line([0 max([XLIM(2) YLIM(2)])], [0 max([XLIM(2) YLIM(2)])]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end




%% ================== RELATIVE DISTANCE FROM NEG VS. POSITIVE CONTROLS
% CHECK MaKE SURE IS CORECT
maxbirds = max(All_birdnum);
maxneur = max(All_neurnum);

% =====================
AllPairMeans = [];
AllBregions = {};
AllBirdnum = [];
for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    for ii=1:maxneur
        
        inds = All_birdnum==i & All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        if strcmp(bregion, 'RA')
            pcol = 'r';
        elseif strcmp(bregion, 'LMAN')
            pcol = 'g';
        end
        
        % ========= COLLECT
        % -- data
        pairtypethis = find(strcmp(PairTypesInOrder, '1  0  0'));
        indstmp = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
        ydat = mean(All_OneMinusRho(indstmp));
               
        % -- neg
        yneg = mean(All_OneMinusRho_NEG(indstmp));
        
        % -- pos
        pairtypethis = find(strcmp(PairTypesInOrder, '1  1  1'));
        indstmp = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
        ypos = mean(All_OneMinusRho(indstmp));
        
        % --- all
        y = (ydat-yneg)/(ypos-yneg);        
        
        % ======== collect each neuron
        AllPairMeans = [AllPairMeans; y];
        AllBregions = [AllBregions; bregion];
        AllBirdnum = [AllBirdnum; i];
        
        % ======== overlay negative distribution
    end
end


% =========== PLOT SCATTER COMPARING TWO CLASSES ACROSS ALL NEURONS
maxbirds = max(AllBirdnum);
for j=1:maxbirds
   lt_figure; hold on;
   xlabel('LMAN -- RA');
   ylabel('position of dat between neg and positive controls')
    Y ={};
   % ---- lman
   inds = AllBirdnum==j & strcmp(AllBregions, 'LMAN');
   Y{1} = AllPairMeans(inds);
   
   % --- ra
   inds = AllBirdnum==j & strcmp(AllBregions, 'RA');
   Y{2} = AllPairMeans(inds);
   
   lt_plot_MultDist(Y, [1 2], 1);
   
    
end


% ========= PLOT OVERALL DISTRIBUTIONS
   lt_figure; hold on;
   title('ALL COMBINED');
   xlabel('LMAN -- RA');
   ylabel('position of dat between neg and positive controls')
    Y ={};
   % ---- lman
   inds = strcmp(AllBregions, 'LMAN');
   Y{1} = AllPairMeans(inds);
   
   % --- ra
   inds = strcmp(AllBregions, 'RA');
   Y{2} = AllPairMeans(inds);
   
   lt_plot_MultDist(Y, [1 2], 1);


