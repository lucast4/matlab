function [OUTSTRUCT, SummaryStruct, Params] = lt_neural_NGRAMS_Compile(dirname, window_prem, Nshuffs)

%% lt 4/24/18 - gets all saved data and extract stats from them and puts into struct

%% ===========
savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' dirname];
Params.dirname = dirname;

%% =========== load metadat
load([savedir '/SummaryStruct.mat']);
load([savedir '/Params.mat']);

%% =========== determine premotor window

% window_prem = [-0.03 0.03]; % relative to syl onset
% Nshuffs = 1; % for negative control

% ------ convert window timing relative to data
motifpredur = Params.regexpr.motifpredur;
windowx = motifpredur+window_prem;
strtype = Params.strtype;
% assert(strcmp(strtype, 'xaa')==1, 'have not coded for other string tyupes yet ...');

%% ========== go thru all birds and extract
All_birdnum =[];
All_neurnum = [];
All_diffsyl_logical = [];
All_OneMinusRho = [];
All_OneMinusRho_NEG = [];
OUTSTRUCT.All_AbsFRdiff_Zrelshuff = [];
OUTSTRUCT.All_AbsFRdiff = [];
OUTSTRUCT.All_AbsFRdiff_NEG = [];
OUTSTRUCT.All_ngrampair_inorder = [];
OUTSTRUCT.All_ngramstring_inorder = {};


numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    tmp = load([savedir '/bird' num2str(i) '.mat']);
    birdstruct = tmp.birdstruct;
    
    % ============================== calculate things across all neurons
    numneurons = length(birdstruct.neuron);
    
    for nn=1:numneurons
        disp([num2str(i) '-' num2str(nn)])
        ngramlist = birdstruct.neuron(nn).ngramlist;
        
        % ----------- all pairwise ngrams
        for j=1:length(ngramlist)
            if isempty(birdstruct.neuron(nn).ngramnum(j).DAT)
                continue
            end
            for jj=j+1:length(ngramlist)
                
                if isempty(birdstruct.neuron(nn).ngramnum(jj).DAT)
                    continue
                end
                
                
                
                % ##################### GET FR MATRIX ACROSS MOTIFS
                % ================ COLLECT ALL TRIALS INTO ONE MATRIX
                if (0) % old version, uses segextract
                    segextract1 = birdstruct.neuron(nn).ngramnum(j).SegmentsExtract;
                    segextract2 = birdstruct.neuron(nn).ngramnum(jj).SegmentsExtract;
                    segextract1 = lt_neural_SmoothFR(segextract1);
                    segextract2 = lt_neural_SmoothFR(segextract2);
                    x = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT.frx;
                    
                    frmat1 = [segextract1.FRsmooth_rate_CommonTrialDur];
                    frmat2 = [segextract2.FRsmooth_rate_CommonTrialDur];
                    if size(frmat1,1)~=size(frmat2,1)
                        indstmp = min([size(frmat1,1), size(frmat2,1)]);
                        frmat1 = frmat1(1:indstmp,:);
                        frmat2 = frmat2(1:indstmp,:);
                    end
                    FRmat = [frmat1 frmat2];
                    TrialInds = {1:length(segextract1), ...
                        length(segextract1)+1:length(segextract1)+length(segextract2)};
                else
                    % new vesrion, directly uses frmat
                    frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
                    frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
                    
                    FRmat = [frmat1 frmat2];
                    n1 = size(frmat1,2);
                    n2 = size(frmat2,2);
                    x = lt_neural_QUICK_XfromFRmat(frmat1);
                    
                    TrialInds = {1:n1, n1+1:n1+n2};
                end
                % --- pare down to just premotor windwo
                FRmat = FRmat(x>=windowx(1) & x<=windowx(2),:);
                
                
                % ############################ TAKE SQUARE ROOT TRANSFORM
                FRmat = sqrt(FRmat);
                
                if (0)
                    % ========== DIAGNOSITCS (HISTOGRAMS)
                    figure; hold on;
                    lt_subplot(2,2,1); hold on;
                    xlabel('noise (spk), mean subtracted at each timebin');
                    tmp = FRmat(:, TrialInds{1}) - mean(FRmat(:, TrialInds{1}),2);
                    lt_plot_histogram(tmp(:));
                    
                    % ---- apply square root transform
                    lt_subplot(2,2,3); hold on;
                    title('sqroot transform each trial');
                    FRmatSq = sqrt(FRmat);
                    tmp = FRmatSq(:, TrialInds{1}) - mean(FRmatSq(:, TrialInds{1}),2);
                    lt_plot_histogram(tmp(:));
                    
                    % ---- apply log transform
                    lt_subplot(2,2,4); hold on;
                    title('log transofrm each trial');
                    FRmatSq = log10(FRmat+0.1*min(FRmat(FRmat(:)>0)));
                    tmp = FRmatSq(:, TrialInds{1}) - mean(FRmatSq(:, TrialInds{1}),2);
                    lt_plot_histogram(tmp(:));
                    
                    % --- variance as function of mean (across time bins)
                    lt_subplot(2,2,2); hold on;
                    xlabel('mean');
                    ylabel('var');
                    mean1 =  mean(FRmat(:, TrialInds{1}),2);
                    var1 = var(FRmat(:, TrialInds{1}),0,2);
                    plot(mean1, var1, 'ok')
                    
                elseif (0)
                    % =============== DIAGNOSTICS (RESIDUAL PLOTS)
                    figure; hold on;
                    lt_subplot(2,2,1); hold on;
                    xlabel('noise (spk), mean subtracted at each timebin');
                    tmp = FRmat(:, TrialInds{1}) - mean(FRmat(:, TrialInds{1}),2);
                    plot(1:size(tmp,1), tmp, '.k')
                    
                    % ---- apply square root transform
                    lt_subplot(2,2,3); hold on;
                    title('sqroot transform each trial');
                    FRmatSq = sqrt(FRmat);
                    tmp = FRmatSq(:, TrialInds{1}) - mean(FRmatSq(:, TrialInds{1}),2);
                    plot(1:size(tmp,1), tmp, '.k')
                    
                    % ---- apply log transform
                    lt_subplot(2,2,4); hold on;
                    title('log transofrm each trial');
                    FRmatSq = log10(FRmat+0.1*min(FRmat(FRmat(:)>0)));
                    tmp = FRmatSq(:, TrialInds{1}) - mean(FRmatSq(:, TrialInds{1}),2);
                    plot(1:size(tmp,1), tmp, '.k')
                    
                    % --- variance as function of mean (across time bins)
                    lt_subplot(2,2,2); hold on;
                    xlabel('mean');
                    ylabel('var');
                    mean1 =  mean(FRmat(:, TrialInds{1}),2);
                    var1 = var(FRmat(:, TrialInds{1}),0,2);
                    plot(mean1, var1, 'ok')
                end
                
                %% ############### RUNNING FR DIFF (relative to shuffles)
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DAT
                % NOTE: have confirmed that distribituions of null data are
                % mostly normal looking. this is becuase of sqrt transfomr
                fr1 = mean(FRmat(:, TrialInds{1}),2);
                fr2 = mean(FRmat(:, TrialInds{2}),2);
                frdiffDAT = mean(abs(fr1-fr2));
                
                OUTSTRUCT.All_AbsFRdiff = ...
                    [OUTSTRUCT.All_AbsFRdiff frdiffDAT];
                %% ########### COMPARE FR DIFF TO SHUFFLE DISTRIBUTION
                FRdiffShuff = [];
                nshufftmp = 100; % 50-100 is optimal, see DIAG below.
                
                nsamps = size(FRmat,2);
                for ss = 1:nshufftmp
                    
                    indshuff = randperm(nsamps);
                    FRmatSHUFF = FRmat(:,indshuff);
                    
                    % ================= calculate correlation
                    fr1 = mean(FRmatSHUFF(:, TrialInds{1}),2);
                    fr2 = mean(FRmatSHUFF(:, TrialInds{2}),2);
                    frdiff = mean(abs(fr1-fr2));
                    
                    FRdiffShuff = [FRdiffShuff frdiff];
                end
                
                % ########################### SAVE THE FIRST SHUFF
                OUTSTRUCT.All_AbsFRdiff_NEG = ...
                    [OUTSTRUCT.All_AbsFRdiff_NEG FRdiffShuff(1)];

                
                
                % ########################## zscore data relative to shuff
                frdiff_Z = (frdiffDAT - mean(FRdiffShuff))./std(FRdiffShuff);
                
                OUTSTRUCT.All_AbsFRdiff_Zrelshuff = ...
                    [OUTSTRUCT.All_AbsFRdiff_Zrelshuff frdiff_Z];
                
                % ======================= DIAG - run this to determine
                % number of shuffles before z-score converges. I think is
                % about 50-100, even for cases of small sample size
                % (n=10-20)
                if (0)
                    Ytmp = [];
                    for nshufftmp = [5 10 15 20 50 75 100 200 500 1000]
                        % ================== COMPARE TO SHUFFLE DISTRIBUTION
                        FRdiffShuff = [];
                        %                 nshufftmp = 200;
                        nsamps = size(FRmat,2);
                        for ss = 1:nshufftmp
                            
                            indshuff = randperm(nsamps);
                            FRmatSHUFF = FRmat(:,indshuff);
                            
                            % ================= calculate correlation
                            fr1 = mean(FRmatSHUFF(:, TrialInds{1}),2);
                            fr2 = mean(FRmatSHUFF(:, TrialInds{2}),2);
                            frdiff = mean(abs(fr1-fr2));
                            
                            FRdiffShuff = [FRdiffShuff frdiff];
                        end
                        
                        % ======================= zscore data relative to shuff
                        y = (frdiffDAT - mean(FRdiffShuff))./std(FRdiffShuff);
                        Ytmp = [Ytmp y];
                    end
                    figure; hold on;
                    subplot(2,2,1); hold on;
                    plot([5 10 15 20 50 75 100 200 500 1000], Ytmp, '-ok');
                    lt_plot_text(100, 1, ['N = ' num2str(n1) ',' num2str(n2)], 'r');
                    subplot(2,2,2); hold on;
                    plot(log10([5 10 15 20 50 75 100 200 500 1000]), Ytmp, '-ok');
                end
                
                
                %% ############################################ RHO
                
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA
                % ================= calculate correlation
                fr1 = mean(FRmat(:, TrialInds{1}),2);
                fr2 = mean(FRmat(:, TrialInds{2}),2);
                rhoDAT = corr(fr1, fr2);
                
                
                if (0) % OLD VERSION, cannot generalize to shuffle analyses.
                    % -------------- get correlation between premotor windows
                    fr1 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT.frmean;
                    x1 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(j).DAT.frx;
                    fr1 = fr1(x1>=windowx(1) & x1<=windowx(2));
                    
                    fr2 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(jj).DAT.frmean;
                    x2 = NGRAMSTRUCT.bird(i).neuron(nn).ngramnum(jj).DAT.frx;
                    fr2 = fr2(x2>=windowx(1) & x2<=windowx(2));
                    rho = corr(fr1, fr2);
                end
                
                
                % ################################## NEGATIVE CONTROL
                RhoShuff = [];
                nsamps = size(FRmat,2);
                for ss = 1:Nshuffs
                    indshuff = randperm(nsamps);
                    
                    FRmatSHUFF = FRmat(:,indshuff);
                    
                    % ================= calculate correlation
                    fr1 = mean(FRmatSHUFF(:, TrialInds{1}),2);
                    fr2 = mean(FRmatSHUFF(:, TrialInds{2}),2);
                    rho = corr(fr1, fr2);
                    
                    RhoShuff = [RhoShuff rho];
                end
                rhoNEG = median(RhoShuff);
                
                % ------------- is this same pair or different pair?
                motif1 = birdstruct.neuron(nn).ngramlist{j};
                motif2 = birdstruct.neuron(nn).ngramlist{jj};
                % -- code: 111 means same at all 3 positions. 010 means
                % diff-same-diff, etc.
                diffsyl_logical = motif1~=motif2;
                
                
                
                %% ################################## OUTPUT
                All_birdnum =[All_birdnum; i];
                All_neurnum = [All_neurnum; nn];
                All_diffsyl_logical = [All_diffsyl_logical; diffsyl_logical];
                All_OneMinusRho = [All_OneMinusRho; 1-rhoDAT];
                All_OneMinusRho_NEG = [All_OneMinusRho_NEG; 1-rhoNEG];
                OUTSTRUCT.All_ngrampair_inorder = [OUTSTRUCT.All_ngrampair_inorder; ...
                    j jj];
                OUTSTRUCT.All_ngramstring_inorder = [OUTSTRUCT.All_ngramstring_inorder; ...
                    {motif1 motif2}];

                
            end
        end
    end
end

%% ==== save in output struct
OUTSTRUCT.All_birdnum =[All_birdnum];
OUTSTRUCT.All_neurnum = [All_neurnum];
OUTSTRUCT.All_diffsyl_logical = [All_diffsyl_logical];
OUTSTRUCT.All_OneMinusRho = [All_OneMinusRho];
OUTSTRUCT.All_OneMinusRho_NEG = [All_OneMinusRho_NEG];


