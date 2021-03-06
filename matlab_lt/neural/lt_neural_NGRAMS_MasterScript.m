%% ================ extract ngrams and fr
close all;

% Params.LearnKeepOnlyBase = 1;
% Params.strtype = 'xaaa';
% Params.Nmin = 6; % num trials minimum
% Params.alignsyl = 3;
% 
% Params.regexpr.motifpredur = 0.2;
% Params.regexpr.motifpostdur = 0.1;
% Params.regexpr.alignOnset = 1;
% Params.regexpr.preAndPostDurRelSameTimept=1;
% Params.regexpr.RemoveIfTooLongGapDur = 1;
Params.LearnKeepOnlyBase = 1;
Params.strtype = 'xaa';
Params.Nmin = 6; % num trials minimum
Params.alignsyl = 2; % which syl to align to...

Params.regexpr.motifpredur = 0.15;
Params.regexpr.motifpostdur = 0.1;
Params.regexpr.alignOnset = 1;
Params.regexpr.preAndPostDurRelSameTimept=1;
Params.regexpr.RemoveIfTooLongGapDur = 1;

saveON =1;

NGRAMSTRUCT = lt_neural_NGRAMS_Extract(SummaryStruct, Params, saveON);



%% #################################### EXTRACTION/PREPROCESSING
%% ================= EXTRACT NGRAMSTRUCT FROM SAVED DATA
close all;
% dirname = 'xaa_30Apr2018_2342'; % 
% dirname = 'xaa_10Mar2019_0109'; % new, from 3/7/2019
dirname = 'xaa_21Mar2019_0146'; % new, from 3/7/2019
window_prem = [-0.025 0.025]; % relative to syl onset
Nshuffs = 1; % for negative control (corr analysis)
doSqrtTransform = 1;
use_dPrime = 0; % if 1, then gets mean dPrime instead of mean abs FR diff
nshufftmp = 2; % for zscoring fr diff.
DoDecode = 0; % to get decoder performances; IMPORTNAT: IF 1, then overwrites
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


%% =========== [PREPROCESS]
% NOTE: currently script, not function, but should run correctly.

lt_neural_NGRAMS_Preprocess;


%% =========== [RE-EXTRACT, EQUALIZING SAMPLE SIZE]
% NOT SURE IF I ACTUALLY DO THIS...
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
DoDecode =0; % IMPORTANT: if this is 1, then uses decode and overwrites FR diff stuff
OUTSTRUCT = lt_neural_NGRAMS_ReSample(OUTSTRUCT, SummaryStruct, Params, ...
    measure_to_recalc, PairTypesToCompare, nshufftmp, DoDecode);

%% ============= save outstruct
savesuffix = '';

fname = ['/bluejay5/lucas/analyses/neural/NGRAMS/' Params.dirname '/OUTSTRUCT_' savesuffix '.mat'];
save(fname, 'OUTSTRUCT');

fname = ['/bluejay5/lucas/analyses/neural/NGRAMS/' Params.dirname '/Params_' savesuffix '.mat'];
save(fname, 'Params');


%% ============= LOAD SAVED STRUCT HERE

load('OUTSTRUCT_.mat');
load('Params_.mat');
load('SummaryStruct.mat');

%% %%%%%%%%%%%%%%%%%%%%%%%% OPTIONAL - DOWNSAMPLING ONE TYPE SO EFFECT SIZE MATCHES OTHER.
% === this script contains ...

lt_neural_NGRAMS_DsampScript;


%% #################################################### [DIAGNOSTICS]
% === coitnained in this script:

lt_neural_NGRAMS_Diagnostic;

%%  ################### DIAGNOSTIC - SEQUENCE OF THINGS I DID TO CHECK [GOOD]

lt_neural_NGRAMS_DiagnScript;


%% ####################################### MAIN PLOTS
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
%     PairtypesToplot = {...
%         '1  1  1', ... % xaxis
%         '0  1  1'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
    PairtypesToplot = {...
        '1  1  1  1', ... % xaxis
        '1  0  0  0'}; % yaxis
end
plotRawGood = 0; % histogram for pos, negative, dat (only works for plottype = absfrdiff_globZ);

removeBadSyls = 1; % i.e. badly labeled...

% ----------- params for one minus rho, specifically
dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
sanitycheckuseneg = 0; % uses negative control data instead of data

% minPairs = 3; % i.e. if any pairtype has fewer pairs that this then will skip
% zscoreMin = [1 0.75]; % [pairnum minz] only keep a neuron if the zscore(mean) for pairnum is greater than minZ
% labelneur = 1;
minPairs = 3; % i.e. if any pairtype has fewer pairs that this then will skip
zscoreMin = [1 0.75]; % [pairnum minz] only keep a neuron if the zscore(mean) for pairnum is greater than minZ
labelneur = 1;

[AllPairs_Means, AllPairs_Birdnum, AllPairs_Bregions] = ...
    lt_neural_NGRAMS_PlotScatter(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg, plotRawGood, usemedian, ...
    removeBadSyls, minPairs, zscoreMin, labelneur);


%% ============ [PLOT] separate by pair type, and average across units
close all;

plottype = 'absfrdiff'; % oneminusrho or absfrdiff
plotON=1; % raw plots? only works for absfrdiff
dosubtractcontrol = 0; % then subtracts negative control before plotting

lt_neural_NGRAMS_PlotByPairtype(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    dosubtractcontrol)


%%  ######################## REGRESSION MODELING

lt_neural_NGRAMS_Regression;


%% ##################################################################
%% ########################### PLOT EXAMPLES
%% ============ for every ngram, figure out how variable gap durations are...


%% ============= PLOT INFORMATION FOR EACH BIRD
% I.E. NGRAMS, SAMPLE SIZES, NEURONS...

pairtype = '1  0  0'; % will list all ngrams that match this pairtype

savedir = ['/bluejay0/bluejay2/lucas/analyses/neural/NGRAMS/' Params.dirname];

lt_neural_NGRAMS_PlotInfo;

%% ============ PLOT RAW NEURAL DATA (multiple trials, for each NGRAM)


close all; 
BirdToPlot = 'gr48bu5';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = [1]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
% motiflist = {'a(b)', 'jbh(h)g'};
% motiflist = {'(d)kcc', 'dk(c)c', '(n)hh', 'c(b)'};
motiflist = {'r(r)d', 'a(r)d'};

% motifpredur = 0.15;
% motifpostdur = 0.15;
motifpredur = 0.1;
motifpostdur = 0.1;
preAndPostDurRelSameTimept = 1;

% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

saveON = 0; 
Nmax = 40;

savedirmain = ['/bluejay0/bluejay2/lucas/analyses/neural/NGRAMS/' Params.dirname '/FIGS/DIAGN_PlotRawNeural/'];

lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
    motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON, ...
    Nmax, savedirmain);



%% ============ [PLOT AND SAVE RAW NEURAL] DOES ALL BIRDS AND NEURONS
close all;

maxbirds = max(OUTSTRUCT.All_birdnum);

for i=1:maxbirds
    
    bname = SummaryStruct.birds(i).birdname;
    
    indsthis = find(OUTSTRUCT.All_birdnum==i);
    
    neurlist = unique(OUTSTRUCT.All_neurnum(indsthis));
    
    ngrams = OUTSTRUCT.All_ngramstring_inorder(indsthis,:);
    ngrams = unique(ngrams(:));
    
    % ===== designate with token syl
    n = Params.alignsyl;
    for j=1:length(ngrams)
        ngrams{j} = [ngrams{j}(1:n-1) '(' ngrams{j}(n) ')' ngrams{j}(n+1:end)];
    end
    
    % ======================== PLOT AND SAVE RAW
    close all;
    BirdToPlot = bname;
    % % ---- give it either
    % A) one neuron and a bunch of motifs or
    % B) bunch of neurons and one motif
    NeurToPlot = neurlist'; % 4 % vector (e.g. [5 7]) - if [] then plots all;
    % motiflist = {'a(b)', 'jbh(h)g'};
    % motiflist = {'(d)kcc', 'dk(c)c', '(n)hh', 'c(b)'};
    motiflist = ngrams';
    
    % motifpredur = 0.15;
    % motifpostdur = 0.15;
    motifpredur = 0.075;
    motifpostdur = 0.075;
    preAndPostDurRelSameTimept = 1;
    
    % --- 1) directed song
    PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both
    
    saveON = 1;
    Nmax = 22;
    savedirmain = ['/bluejay0/bluejay2/lucas/analyses/neural/NGRAMS/' Params.dirname '/FIGS/DIAGN_PlotRawNeural/'];

    lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
        motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON, ...
        Nmax, savedirmain);
    
end


%% =============== [RAW NEURAL] After saving (above) then plot
birdtoplot = 'pu69wh78';
neurtoplot = 11;
motiftoplot ='n(a)a';

% close all;

% dname = dir('*');
% dname = dir(['wh44wh39-neur*' '-' motiftoplot]);
dname = dir([birdtoplot '-neur' num2str(neurtoplot) '-' motiftoplot]);

for i=1:length(dname)
    disp(i);
    if dname(i).isdir
        try
        openfig([dname(i).name '/figsall.fig']);
         
        set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
%         figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        ylim([-200 200]);
        pause; 
        catch err
            
        end

    end
    
    if mod(i,10)==0
        close all; 
    end
    if mod(i,11)==0
        close all; 
    end
end

%% =============== [RAW NEURAL] After saving (above) then plot
birdtoplot = 'bk7';
neurtoplot = 10;
motiftoplot ='g(b)b';

openfig([birdtoplot '-neur' num2str(neurtoplot) '-' motiftoplot '/figsall.fig']);
% 
% for i=1:length(dname)
%     disp(i);
%     if dname(i).isdir
%         try
%         openfig([dname(i).name '/figsall.fig']);
%          
%         set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% %         figure('units', 'normalized', 'outerposition', [0 0 1 1]);
%         ylim([-200 200]);
%         pause; 
%         catch err
%             
%         end
% 
%     end
%     
%     if mod(i,10)==0
%         close all; 
%     end
%     if mod(i,11)==0
%         close all; 
%     end
% end
% 
% 

%% ============= [RAW NEURAL] one plot for each neuron (i.e. to check spike sorting
% GO TO SAVE FOLDER FROM ABOVE.
cc =1 ;
for i=1:length(SummaryStruct.birds)
    for ii=1:length(SummaryStruct.birds(i).neurons)
        disp(['bird' num2str(i) ' - neur' num2str(ii)]);
        bname = SummaryStruct.birds(i).birdname;
        
        fname = [bname '-neur' num2str(ii) '*'];
        tmp = dir(fname);
        
        if isempty(tmp)
            continue
        end
        
        openfig([tmp(1).name '/figsall.fig']);
        
        set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
        %         figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        ylim([-200 200]);
        pause;
        cc=cc+1;
        if mod(cc,10)==0
            close all;
        end
    end
end

%% #############################################################
%% ########################## EXAMPLES, V2 - good
%% ============= [SHOW METADAT] FOR EACH BIRD AND NEURON, SHOW DIRECTORY
% for each directory, displays channels to get.
longversion=1;
lt_neural_SummaryStruct_ShowDirChans(SummaryStruct, longversion);


%% ============= PLOT EXAMPLE FR TRACES
% only plots 2. below plots more (but with less detail)
if (0)
    close all;
birdtoplot = 'wh44wh39';
neurtoplot = 71;
ngramstoplot = {}; % leave empty to plot random one
pairtypetoplot = '1  1  1';

lt_neural_NGRAMS_PlotEgAll(OUTSTRUCT, SummaryStruct, Params, ...
    birdtoplot, neurtoplot, ngramstoplot, pairtypetoplot);

end

%% ============== PLOT ALL EXAMPLES FOR A GIVEN NEURON IN ONE PLOT

close all;
birdtoplot = 'gr48bu5';
neurtoplot = 18;
pairtypes = {'1  0  0', '1  1  1'};
plotsqrt = 0; % if 1, then mean of sqrt firing rate
lt_neural_NGRAMS_PlotEgPair(OUTSTRUCT, SummaryStruct, Params, ...
    birdtoplot, neurtoplot, pairtypes, plotsqrt);



%% =========== EXAMPLE - plot motif pair strings...
% close all;

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
neurtoplot = 1;

lt_neural_NGRAMS_PlotMotStr(OUTSTRUCT, SummaryStruct, plottype, ...
    PairtypesToplot, birdtoplot, neurtoplot);



%% ##################################### [OLDER PLOTS]
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




%% ################################ TIMECOURSE ANALYSIS
%% EXTRACT/COMPUTE TIMECOURSES FROM SAVED RAW DAT

close all;
measure_to_recalc = 'absfrdiff';
if strcmp(Params.strtype, 'xaa')
%     PairTypesToCompare = {...
%         '1  1  1', ... % xaxis
%         '1  0  0', ...
%         '0  0  1'}; % yaxis
    PairTypesToCompare = {...,
        '0  0  1', ...
        '1  0  0', ...
        '1  0  1', ...
        '0  1  0', ...
        '0  1  1', ...
        '1  1  0', ...
        '1  1  1'}; % in order to be plotted
elseif strcmp(Params.strtype, 'xaaa')
    PairTypesToCompare = {...
        '1  1  1  1', ... % xaxis
        '1  0  0  0'}; % yaxis
end
nshufftmp = 2;
DoDecode =1; % IMPORTANT: if this is 1, then uses decode and overwrites FR diff stuff
tcoursestyle = 'time';
% time: diff in FR (1ms bins)
% binned: overlapping windows [IN PROGRESS];
OUTSTRUCT = lt_neural_NGRAMS_Timecourse(OUTSTRUCT, SummaryStruct, Params, ...
    PairTypesToCompare, nshufftmp, DoDecode, tcoursestyle);


% =================== SAVE
savesuffix = 'Tcourse';
fname = ['/bluejay5/lucas/analyses/neural/NGRAMS/' Params.dirname '/OUTSTRUCT_' savesuffix '.mat'];
save(fname, 'OUTSTRUCT');
fname = ['/bluejay5/lucas/analyses/neural/NGRAMS/' Params.dirname '/Params_' savesuffix '.mat'];
save(fname, 'Params');



%% ========================== COMPUTE RUNNING CONTEXTUAL MODULATION INDEX

close all;
if strcmp(Params.strtype, 'xaa')
    PairTypesToCompare = {...
        '1  1  1', ... % xaxis
        '1  0  0'}; % yaxis
    PairTypesToCompare = {...
        '1  1  1', ... % xaxis
        '0  0  1'}; % yaxis
    PairTypesToCompare = {...
        '1  1  1', ... % xaxis
        '0  1  1'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
    PairTypesToCompare = {...
        '1  1  1  1', ... % xaxis
        '1  0  0  0'}; % yaxis
end
tcoursestyle = 'time';
useGlobalNeg = 0;
numrawplotsperbird = 0;
lt_neural_NGRAMS_Timecourse_Calc(OUTSTRUCT, SummaryStruct, Params, ...
    PairTypesToCompare, tcoursestyle, useGlobalNeg, numrawplotsperbird)




%% #############################################################
%% ########################### FF AND PITCH CORRELATIONS ANALYSES



















