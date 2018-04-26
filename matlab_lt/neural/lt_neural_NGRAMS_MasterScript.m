%% ================ extract ngrams and fr
close all;

Params.LearnKeepOnlyBase = 1;
Params.strtype = 'xaa';
Params.Nmin = 6; % num trials minimum
Params.alignsyl = 2;

Params.regexpr.motifpredur = 0.1;
Params.regexpr.motifpostdur = 0.1;
Params.regexpr.alignOnset = 1;
Params.regexpr.preAndPostDurRelSameTimept=1;
Params.regexpr.RemoveIfTooLongGapDur = 1;

saveON =1;

NGRAMSTRUCT = lt_neural_NGRAMS_Extract(SummaryStruct, Params, saveON);


%% ================= EXTRACT NGRAMSTRUCT FROM SAVED DATA 
close all;
dirname = 'xaa_26Apr2018_0938';
window_prem = [-0.025 0.025]; % relative to syl onset
Nshuffs = 1; % for negative control (corr analysis)

[OUTSTRUCT, SummaryStruct, Params] = lt_neural_NGRAMS_Compile(dirname, ...
    window_prem, Nshuffs);

% ========================== FOR COMPATIBILITY WITH OLD CODE, EXTRACT ALL
% FIELDS
fnamesthis = fieldnames(OUTSTRUCT);
for j=1:length(fnamesthis)
   eval([fnamesthis{j} ' = OUTSTRUCT.' fnamesthis{j} ';']); 
end

%% ================= get all pairwise distances during premotor window
if (0) % OLD VERSION --- this works with NGRAMSTRUCT. new version does not since 
    % filesize too large.
    
    % if ever want to use this version, need to run the script in here.
    % this extracts summary arrays.
lt_neural_NGRAMS_GetDatOldVersion;

end
%% ================== convert motifpair types to groups

PairTypesInOrder = {...,
    '0  0  1', ...
    '1  0  0', ...
    '1  0  1', ...
    '0  1  0', ...
    '0  1  1', ...
    '1  1  0', ...
    '1  1  1'}; % in order to be plotted

All_diffsyl_string = num2str(double(All_diffsyl_logical));
All_diffsyl_string = mat2cell(All_diffsyl_string, ones(size(All_diffsyl_logical,1),1));

[~, All_diffsyl_PairType] = ismember(All_diffsyl_string, PairTypesInOrder);
assert(all(strcmp(All_diffsyl_string, PairTypesInOrder(All_diffsyl_PairType)')), 'asdfas');

% ================= PUT INTO STRUCT
OUTSTRUCT.PairTypesInOrder = PairTypesInOrder;
OUTSTRUCT.All_diffsyl_PairType = All_diffsyl_PairType;

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

plottype = 'absfrdiff'; % oneminusrho or absfrdiff
plotON=0; % only works for absfrdiff
PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis

% ----------- params for one minus rho, specifically
dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
sanitycheckuseneg = 0; % uses negative control data instead of data

lt_neural_NGRAMS_PlotScatter(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg);



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


