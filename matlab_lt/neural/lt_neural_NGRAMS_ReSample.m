function [OUTSTRUCT, Params] = lt_neural_NGRAMS_ReSample(OUTSTRUCT, SummaryStruct, Params, ...
    measure_to_recalc, PairTypesToCompare, nshufftmp, DoDecode)
% nshufftmp = 10; % 50-100 is optimal, see DIAG below.
% NOTE: usually 100, but only matters for data in which
% want to zscore vs. self. for global z this just takes the
% first shuffle anyways.

% PairTypesToCompare = {'1  0  0', '1  1  1'};
% measure_to_recalc = 'absfrdiff'; % currently only 'absfrdiff' works

%% lt 4/30/18 - given 2 pairtypes, resample higher sample pairtype to equalize sample size distributions

% Given two pairtypes, downsamples one pairtype such that the mean sample
% size is the same between pairtypes.

% Also checks whether negative shuffle control magnitude is a function of
% sample size

% Also confirms that negative shuffle control does not change.

% Also plots data as a function of resampled N

%%

% ------ convert window timing relative to data
motifpredur = Params.regexpr.motifpredur;
windowx = motifpredur+Params.window_prem;
strtype = Params.strtype;
% assert(strcmp(strtype, 'xaa')==1, 'have not coded for other string tyupes yet ...');

%% ===========
savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' Params.dirname];


%% ========= note down in params

Params.DidResample=1;

%% ========== go thru all birds and extract

% ============================= FIRST, move old version
if strcmp(measure_to_recalc, 'absfrdiff')
    OUTSTRUCT.All_AbsFRdiff_Zrelshuff_ORIG = OUTSTRUCT.All_AbsFRdiff_Zrelshuff;
    OUTSTRUCT.All_AbsFRdiff_ORIG = OUTSTRUCT.All_AbsFRdiff;
    OUTSTRUCT.All_AbsFRdiff_NEG_ORIG = OUTSTRUCT.All_AbsFRdiff_NEG;
    
    OUTSTRUCT.All_AbsFRdiff_Zrelshuff = [];
    OUTSTRUCT.All_AbsFRdiff = [];
    OUTSTRUCT.All_AbsFRdiff_NEG = [];
    
    OUTSTRUCT.All_N_ORIG = [];
    OUTSTRUCT.All_N = [];
end


% ===================== SECOND, RUN
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    tmp = load([savedir '/bird' num2str(i) '.mat']);
    birdstruct = tmp.birdstruct;
    birdname = SummaryStruct.birds(i).birdname;
    
    % ============================== calculate things across all neurons
    numneurons = length(birdstruct.neuron);
    
    for nn=1:numneurons
        disp(' ############################################################### ');
        disp([birdname ', ' num2str(i) '-' num2str(nn)]);
        ngramlist = birdstruct.neuron(nn).ngramlist;
        
        MotifPairs = {};
        SampleSizes =[];
        PairType = [];
        
        % ----------- all pairwise ngrams
        for j=1:length(ngramlist)
            if isempty(birdstruct.neuron(nn).ngramnum(j).DAT)
                continue
            end
            for jj=j+1:length(ngramlist)
                
                if isempty(birdstruct.neuron(nn).ngramnum(jj).DAT)
                    continue
                end
                
                % ================ skip this if have already removed
                % because is a bad syl
                % ---- pairtypes
                indtmp = find(OUTSTRUCT.All_birdnum == i & OUTSTRUCT.All_neurnum==nn & ...
                    all(OUTSTRUCT.All_ngrampair_inorder == [j jj],2));
                if isempty(indtmp)
                    % then this has been removed
                    continue
                end
                
                assert(length(indtmp)<2, 'OUTSTRUCT and other dat arent matched');
                pairtype = OUTSTRUCT.All_diffsyl_PairType(indtmp);
                
                
                
                % =================== 1) FOR EACH PAIR TY
                % ---- syl 1
                motif1 = birdstruct.neuron(nn).ngramnum(j).regexprstr;
                motif2 = birdstruct.neuron(nn).ngramnum(jj).regexprstr;
                
                % ---- sample size
                N1 = size(birdstruct.neuron(nn).ngramnum(j).DAT.frmat, 2);
                N2 = size(birdstruct.neuron(nn).ngramnum(jj).DAT.frmat, 2);
                
                % ==================== OUTPUT
                % ------- this neuron
                MotifPairs = [MotifPairs; {motif1, motif2}];
                SampleSizes = [SampleSizes; [N1 N2]];
                PairType = [PairType; pairtype];
                
                % ------ all data
                %                 All_MotifPairs = [All_MotifPairs; {motif1, motif2}];
                %                 All_SampleSizes = [All_SampleSizes; [N1 N2]];
                %                 All_PairType = [All_PairType; pairtype];
                %                 All_Birdnum = [All_Birdnum; i];
                %                 All_Neurnum = [All_Neurnum; nn];
                
                % ============================== OUTPUT SOME THINGS
                
                
                %% ======= [SANITY CHECK] frdiff (DAT and SHUFF) as function of sample size (using resampling)
                if (0)
                    if N1>100 & N2>100
                        
                        % ============================= see what shuff fr diff
                        % is as a function of sample size
                        frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
                        frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
                        
                        FRmat = [frmat1 frmat2];
                        x = lt_neural_QUICK_XfromFRmat(frmat1);
                        TrialInds = {1:N1, N1+1:N1+N2};
                        
                        % -------------- pare down to just premotor windwo
                        FRmat = FRmat(x>=windowx(1) & x<=windowx(2),:);
                        % ------------------- TAKE SQUARE ROOT TRANSFORM
                        if Params.doSqrtTransform==1
                            FRmat = sqrt(FRmat);
                        end
                        
                        % ############### RUNNING FR DIFF (relative to shuffles)
                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DAT
                        nlist = [6 8 10 14 18 20 25 30 40 50 75 100 150 200];
                        FRDatmat = nan(length(nlist), length(nlist)); % holds actual diff
                        FRShuffmat = nan(length(nlist), length(nlist)); % holds shuffle
                        for k = 1:length(nlist)
                            z = nlist(k);
                            for kk = 1:length(nlist)
                                zz = nlist(kk);
                                
                                % --- get random inds and do subssample
                                indperm1 = randperm(N1, z);
                                indperm2 = randperm(N2, zz);
                                
                                fr1 = mean(FRmat(:, TrialInds{1}(indperm1)),2);
                                fr2 = mean(FRmat(:, TrialInds{2}(indperm2)),2);
                                frdiffDAT = mean(abs(fr1-fr2));
                                
                                FRDatmat(k, kk) = frdiffDAT;
                                
                                % ----- get shuffle
                                FRmat_shufftmp = FRmat(randperm(size(FRmat,2)));
                                fr1 = mean(FRmat_shufftmp(:, TrialInds{1}(indperm1)),2);
                                fr2 = mean(FRmat_shufftmp(:, TrialInds{2}(indperm2)),2);
                                frdiffDAT = mean(abs(fr1-fr2));
                                
                                FRShuffmat(k, kk) = frdiffDAT;
                            end
                        end
                        
                        % ============= plot
                        lt_figure; hold on;
                        % --- DAT
                        lt_subplot(2,2,1); hold on;
                        title('[DAT] diff N1, diff colors');
                        xlabel('N2');
                        ylabel('frdiff');
                        pcols = lt_make_plot_colors(length(nlist), 1, [1 0 0]);
                        for k=1:length(nlist)
                            frall = FRDatmat(k,:);
                            plot(nlist, frall, '-o', 'Color', pcols{k});
                        end
                        
                        lt_subplot(2,2,2); hold on;
                        title('[DAT] mean N between motifs');
                        xlabel('mean N');
                        ylabel('fr diff');
                        
                        x = [];
                        y = [];
                        for k=1:length(nlist)
                            for kk=1:length(nlist)
                                x = [x mean([nlist(k) nlist(kk)])];
                                y = [y FRDatmat(k, kk)];
                            end
                        end
                        plot(x, y, 'ok');
                        Xdat = x;
                        Ydat = y;
                        
                        % --- SHUFFLE
                        lt_subplot(2,2,3); hold on;
                        title('[SHUFF] diff N1, diff colors');
                        xlabel('N2');
                        ylabel('frdiff');
                        pcols = lt_make_plot_colors(length(nlist), 1, [1 0 0]);
                        for k=1:length(nlist)
                            frall = FRShuffmat(k,:);
                            plot(nlist, frall, '-o', 'Color', pcols{k});
                        end
                        
                        lt_subplot(2,2,4); hold on;
                        title('[SHUFF] mean N between motifs');
                        xlabel('mean N');
                        ylabel('fr diff');
                        
                        x = [];
                        y = [];
                        for k=1:length(nlist)
                            for kk=1:length(nlist)
                                x = [x mean([nlist(k) nlist(kk)])];
                                y = [y FRShuffmat(k, kk)];
                            end
                        end
                        plot(x, y, 'ok');
                        
                        Xshuff = x;
                        Yshuff = y;
                        
                        % ----
                        lt_figure; hold on
                        title('each point different resampled size');
                        xlabel('dat');
                        ylabel('shuff');
                        lt_plot_45degScatter(Ydat, Yshuff);
                    end
                end
            end
        end
        
        % ==========================
        if isempty(PairType)
            continue
        end
        
        %% ======= [FOR THIS NEURON] figure out relative distriubtions of sample sizes
        Indpaircomp = find(ismember(OUTSTRUCT.PairTypesInOrder, PairTypesToCompare));
        
        % ======= Get distribitiosn of sample sizes
        indsType1 = PairType == Indpaircomp(1);
        indsType2 = PairType == Indpaircomp(2);
        
        Nall = {};
        Nall{1} = SampleSizes(indsType1, :);
        Nall{2} = SampleSizes(indsType2, :);
        
        % ====== which type has larger sample size?
        if mean(Nall{2}(:)) > mean(Nall{1}(:))
            ind_larger = 2;
        elseif mean(Nall{1}(:)) > mean(Nall{2}(:))
            ind_larger = 1;
        else
            ind_larger = []; % then they are already identical...
        end
        assert(~isempty(ind_larger), 'why exaclty same sample size?');
        
        
        % =========== save: which pairtype to downsample and by how much
        % --- which pairtype
        pairtypeToDwnsmp = Indpaircomp(ind_larger);
        
        % --- compare how much in excess of the Nmin
        nexcess_1 = mean(Nall{1}(:)) - Params.Nmin;
        nexcess_2 = mean(Nall{2}(:)) - Params.Nmin;
        downfactor = [];
        if ind_larger ==1
            downfactor = nexcess_2/nexcess_1;
        elseif ind_larger ==2
            downfactor = nexcess_1/nexcess_2;
        end
        
        
        %% ========== re calculate fr diff, this time after downsampling
        % ----------- all pairwise ngrams
        for j=1:length(ngramlist)
            if isempty(birdstruct.neuron(nn).ngramnum(j).DAT)
                continue
            end
            for jj=j+1:length(ngramlist)
                
                if isempty(birdstruct.neuron(nn).ngramnum(jj).DAT)
                    continue
                end
                
                % ===== skip this if has been thrown out because is b ad
                % syl
                % --- 1) is the current pairtype the one that should be
                % downsampled?
                indtmp = find(OUTSTRUCT.All_birdnum == i & OUTSTRUCT.All_neurnum==nn & ...
                    all(OUTSTRUCT.All_ngrampair_inorder == [j jj],2));
                if isempty(indtmp)
                    continue
                end
                
                assert(length(indtmp)<2, 'OUTSTRUCT and other dat arent matched');
                pairtype = OUTSTRUCT.All_diffsyl_PairType(indtmp);
                pairtype_logical = OUTSTRUCT.All_diffsyl_logical(indtmp, :);
                
                
                % ============ GET FRMAT across 2 motifs
                frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
                frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
                N1 = size(frmat1,2);
                N2 = size(frmat2, 2);
                assert(N1>=Params.Nmin, 'asfsd');
                assert(N2>=Params.Nmin, 'asfsd');
                
                
                % ================== OUTPUT SOME THINGS
                OUTSTRUCT.All_N_ORIG = [OUTSTRUCT.All_N_ORIG; [N1 N2]];
                
                
                
                %% =========== get fr differences
                if pairtype == pairtypeToDwnsmp
                    % then do downsample
                    [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds, ...
                    DecodeDAT, DecodeShuff] = lt_neural_NGRAMS_QUICK_FRdiff(...
                        frmat1, frmat2, downfactor, windowx, Params.Nmin, nshufftmp, ...
                        Params, DoDecode);
                else
                    % don't do downsample
                    [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds, ...
                    DecodeDAT, DecodeShuff] = lt_neural_NGRAMS_QUICK_FRdiff(...
                        frmat1, frmat2, 1, windowx, Params.Nmin, nshufftmp, Params, ...
                        DoDecode);
                end
                
                OUTSTRUCT.All_N = [OUTSTRUCT.All_N; Nboth'];
                
                
                OUTSTRUCT.All_AbsFRdiff = ...
                    [OUTSTRUCT.All_AbsFRdiff; FRdiffDAT];
                OUTSTRUCT.All_AbsFRdiff_NEG = ...
                    [OUTSTRUCT.All_AbsFRdiff_NEG; FRdiffShuff(1)];
                OUTSTRUCT.All_AbsFRdiff_Zrelshuff = ...
                    [OUTSTRUCT.All_AbsFRdiff_Zrelshuff; FRdiff_Z];
                
                
                
                
                %%
                if (0) % OLDER VERSION...
                    % ============ DOWNSAMPLE
                    
                    if pairtype == pairtypeToDwnsmp
                        % then downsample
                        N1_new = round(downfactor*(N1 - Params.Nmin)) + Params.Nmin; % downfactor*[excess from Nmin] + [Nmin]
                        N2_new = round(downfactor*(N2 - Params.Nmin)) + Params.Nmin; % downfactor*[excess from Nmin] + [Nmin]
                        
                        % ================= replace things
                        frmat1 = frmat1(:, randperm(N1, N1_new));
                        frmat2 = frmat2(:, randperm(N2, N2_new));
                        N1 = N1_new;
                        N2 = N2_new;
                    end
                    
                    % =========== OUTPUT STUFF
                    OUTSTRUCT.All_N = [OUTSTRUCT.All_N; [N1 N2]];
                    
                    
                    % =========== Concatenate into FRmat
                    FRmat = [frmat1 frmat2];
                    x = lt_neural_QUICK_XfromFRmat(frmat1);
                    TrialInds = {1:N1, N1+1:N1+N2};
                    
                    % -------------- pare down to just premotor windwo
                    FRmat = FRmat(x>=windowx(1) & x<=windowx(2),:);
                    % ------------------- TAKE SQUARE ROOT TRANSFORM
                    if Params.doSqrtTransform==1
                        FRmat = sqrt(FRmat);
                    end
                    
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DAT
                    fr1 = mean(FRmat(:, TrialInds{1}),2);
                    fr2 = mean(FRmat(:, TrialInds{2}),2);
                    frdiffDAT = mean(abs(fr1-fr2));
                    
                    OUTSTRUCT.All_AbsFRdiff = ...
                        [OUTSTRUCT.All_AbsFRdiff; frdiffDAT];
                    
                    % ############################ CONTROLS
                    FRdiffShuff = [];
                    
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
                        [OUTSTRUCT.All_AbsFRdiff_NEG; FRdiffShuff(1)];
                    
                    
                    % ########################## zscore data relative to shuff
                    frdiff_Z = (frdiffDAT - mean(FRdiffShuff))./std(FRdiffShuff);
                    
                    OUTSTRUCT.All_AbsFRdiff_Zrelshuff = ...
                        [OUTSTRUCT.All_AbsFRdiff_Zrelshuff; frdiff_Z];
                    
                end
            %% ======================= SANITY CHECK
            motif1 = birdstruct.neuron(nn).ngramlist{j};
            motif2 = birdstruct.neuron(nn).ngramlist{jj};
            % -- code: 111 means same at all 3 positions. 010 means
            % diff-same-diff, etc.
            diffsyl_logical = motif1~=motif2;
            assert(all(diffsyl_logical == pairtype_logical), 'not matched ...');
            end
            
            
            
        end
    end
end
end



















