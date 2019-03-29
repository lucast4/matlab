%% ######################################################### 
%% ################### CHECK NEURAL RAW FOR MOVEMENT ARTIFACT
% THIS IS GENERALLY APPLICABLE, NOT JUST FOR NGRAMS
% TO GO THRU ALL DIRECTORIES BY HAND, THIS TELLS YOU WHERE TO GO:

% ============= [SHOW METADAT] FOR EACH BIRD AND NEURON, SHOW DIRECTORY
% for each directory, displays channels to get.
longversion=0;
lt_neural_SummaryStruct_ShowDirChans(SummaryStruct, longversion);


%% runt his to quickly go to each directory by hand
tmp = '/bluejay0/bluejay2/lucas/birds/gr48bu5/NEURAL/012819_RALMANLearn6  ---  8  15  20  23'
indtmp = strfind(tmp, '---');
cd(tmp(1:indtmp-3));
eval(['!ls Batch*'])

%% ==== OPEN BATCH FILE, CHECK A FEW SONGS. NOTE DOWN BAD THINGS IN 
lt_neural_QUICK_RemoveBadChans;
lt_neural_QUICK_RemoveBadSyl;
lt_neural_QUICK_RemoveTrials;


%% ######################################################### 
%% ################### CHECK THAT SPIKE SORTING IS REASONABLE.


%% [RAW NEURAL] one plot for each neuron (i.e. to check spike sorting
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


%% ######################################################### 
%% ################### TO CHECK RAW NEURAL - DID THIS ONLY FOR NGRAMS THAT WERE PART OF CONVERGENT BRANCH POINT
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


%% ============= LIST ALL NGRAMS 
%% ============= PLOT INFORMATION FOR EACH BIRD
% I.E. NGRAMS, SAMPLE SIZES, NEURONS...

pairtype = '1  0  0'; % will list all ngrams that match this pairtype


savedir = ['/bluejay0/bluejay2/lucas/analyses/neural/NGRAMS/' Params.dirname];

lt_neural_NGRAMS_PlotInfo;


%% FOR A GIVEN BIRD, PLOT ALL NEURONS FOR A GIVEN NGRAM...
% TO FIGUR EOUT WHICH NGRAMS, RUN ABOVE FIRST

% =============== [RAW NEURAL] After saving (above) then plot
birdtoplot = 7;
motiftoplot ='a(j)j';
dname = dir('*');
dname = dir(['wh44wh39-neur*' '-' motiftoplot]);

% ========= RUN
close all;
for i=1:2:length(dname)
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

%% ######################################################### 
%% ################### TO CHECK PSTH, ALL , FOR A GIVEN NEURON.
%% ================ 1) PLOT ALL NEURONS
% =========== [SCATTER PLOT] for each bird, plot scatter of two pairtypes

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

minPairs = 3; % i.e. if any pairtype has fewer pairs that this then will skip
zscoreMin = [1 0.75]; % [pairnum minz] only keep a neuron if the zscore(mean) for pairnum is greater than minZ

labelneur = 1;

[AllPairs_Means, AllPairs_Birdnum, AllPairs_Bregions] = ...
    lt_neural_NGRAMS_PlotScatter(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg, plotRawGood, usemedian, ...
    removeBadSyls, minPairs, zscoreMin, labelneur);

%% ===================== 2) PLOT ALL EXAMPLES FOR A GIVEN NEURON IN ONE PLOT

% close all;
birdtoplot = 'bk7';
neurtoplot = 10;
pairtypes = {'1  0  0', '1  1  1'};
plotsqrt = 1; % if 1, then mean of sqrt firing rate
lt_neural_NGRAMS_PlotEgPair(OUTSTRUCT, SummaryStruct, Params, ...
    birdtoplot, neurtoplot, pairtypes, plotsqrt);

%% ===================== 3) if required, check raw data:
birdtoplot = 'bk7';
neurtoplot = 10;
motiftoplot ='g(b)b';

openfig([birdtoplot '-neur' num2str(neurtoplot) '-' motiftoplot '/figsall.fig']);

