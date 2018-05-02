function lt_neural_NGRAMS_DIAGNOSTIC(OUTSTRUCT, SummaryStruct, Params, ...
    dispNgramStrings, plotRawSampSize, compareNtoFRdist)
%% lt 4/26/18 - sanity checks

%% % ============= 1) SAMPLE SIZE (TRIALS) DIFFERENCES DEPENDING ON PAIRTYPE?
% ----------------- FOR EACH UNIT PLOT DISTRIBUTIONS OF SAMPLE TYPES.

savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' Params.dirname];



%% ========== go thru all birds and extract
All_MotifPairs = {};
All_SampleSizes =[];
All_PairType = [];
All_Birdnum = [];
All_Neurnum = [];

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
                
                % =================== 1) FOR EACH PAIR TY
                % ---- syl 1
                motif1 = birdstruct.neuron(nn).ngramnum(j).regexprstr;
                motif2 = birdstruct.neuron(nn).ngramnum(jj).regexprstr;
                
                % ---- sample size
                N1 = size(birdstruct.neuron(nn).ngramnum(j).DAT.frmat, 2);
                N2 = size(birdstruct.neuron(nn).ngramnum(jj).DAT.frmat, 2);
                
                % ---- pairtypes
                indtmp = find(OUTSTRUCT.All_birdnum == i & OUTSTRUCT.All_neurnum==nn & ...
                    all(OUTSTRUCT.All_ngrampair_inorder == [j jj],2));
                assert(length(indtmp)==1, 'OUTSTRUCT and other dat arent matched');
                pairtype = OUTSTRUCT.All_diffsyl_PairType(indtmp);
                
                
                % ==================== OUTPUT
                % ------- this neuron
                MotifPairs = [MotifPairs; {motif1, motif2}];
                SampleSizes = [SampleSizes; [N1 N2]];
                PairType = [PairType; pairtype];
                
                % ------ all data
                All_MotifPairs = [All_MotifPairs; {motif1, motif2}];
                All_SampleSizes = [All_SampleSizes; [N1 N2]];
                All_PairType = [All_PairType; pairtype];
                All_Birdnum = [All_Birdnum; i];
                All_Neurnum = [All_Neurnum; nn];
                
            end
        end
        
        %% ====== List ngram pairs, grouped by pairtype?
        if dispNgramStrings==1
            if rand<0.1
                % ----- display pairtypes
                maxpairtype = max(PairType);
                for p = 1:maxpairtype
                    
                    % ==== display all pairs of this type
                    indtmp = PairType == p;
                    ptypethis = OUTSTRUCT.PairTypesInOrder{p};
                    
                    ngrampairs = MotifPairs(indtmp,:);
                    
                    % --- if ngrampairs too large N, then get random subsample
                    if size(ngrampairs,1)>200
                        indperm = randperm(size(ngrampairs,1), 200);
                        ngrampairs = ngrampairs(indperm,:);
                    end
                    
                    disp('============================');
                    disp(['Pair type ' num2str(p) ' [ ' ptypethis ' ] N = ' num2str(sum(indtmp))]);
                    disp(ngrampairs');
                    
                end
            end
        end
        
    end
end


%% ====== plot sample size distrubutions across all pairtypes?
PairTypesToCompare = {'1  0  0', '1  1  1'};
indpaircomp = find(ismember(OUTSTRUCT.PairTypesInOrder, PairTypesToCompare));
assert(length(indpaircomp)==2,'asdf');
pcols_pairs = lt_make_plot_colors(length(indpaircomp), 0,0);

if plotRawSampSize==1
    
    figcount=1;
    subplotrows=6;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    maxbirds = max(All_Birdnum);
    maxneur = max(All_Neurnum);
    
    NComp_birdnum = [];
    NComp_neurnum= [];
    NComp_pval= [];
    NComp_N= [];
    
    
    for i=1:maxbirds
        birdname = SummaryStruct.birds(i).birdname;
        for ii=1:maxneur
            
            inds = All_Birdnum==i & All_Neurnum==ii;
            
            if ~any(inds)
                continue
            end
            
            pairtypes = All_PairType(inds);
            sampsize = All_SampleSizes(inds,:);
            
            
            % ===== plot paired lines for sample sizes
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname ',' num2str(i) '-' num2str(ii)]);
            ylabel('n trials');
            xlabel('pairtype (N capped at 200)');
            
            for k =1:max(pairtypes)
                
                indtmp = pairtypes == k;
                if ~any(indtmp)
                    continue
                end
                
                x = [k-0.3 k+0.3];
                y = sampsize(indtmp, :);
                
                
                % --- if too large, pare down sample size
                if size(y,1) > 200
                    indperm = randperm(size(y,1), 200);
                    
                    y = y(indperm,:);
                end
                
                % ------- PLOT
                plot(x, y', '--', 'Color', [0.7 0.7 0.7]);
                
                ymean = mean(y,1);
                ystd = std(y, 0, 1);
                
                lt_plot(x, ymean, {'Errors', ystd});
                
                
                % ------- plot mean across ngrams (i.e. each pair gets one
                % N)
                ymean_comb = mean(mean(y,2),1);
                ystd_comb = std(mean(y,2), 0,1);
                lt_plot(k, ymean_comb, {'Errors', ystd_comb, 'Color', [0.2 0.8 0.2], 'LineStyle', '-'});
                
                
            end
            
            set(gca, 'XTick', 1:max(pairtypes));
            xlim([0 max(pairtypes)+1]);
            lt_plot_zeroline;
            
           %% ======= do statistical comparison between 2 pairtypes
            Y = {}; % to store sample sizes
            
            % --- go through both pairs
            for pnum = [1 2]
                
                indtmp = All_Birdnum==i & All_Neurnum==ii & ...
                    All_PairType==indpaircomp(pnum);
                
                Y{pnum} = mean(All_SampleSizes(indtmp,:), 2); % mean of 2 ngrams
            end
            
            p = ranksum(Y{1}, Y{2});
            lt_plot_pvalue(p, ...
                ['ranksum (' num2str(indpaircomp(1)) 'vs' num2str(indpaircomp(2)) ')'], 1);
            
            
            % ================= OUTPUT
            NComp_birdnum = [NComp_birdnum; i];
            NComp_neurnum = [NComp_neurnum; ii];
            NComp_pval = [NComp_pval; p];
            NComp_N = [NComp_N; cellfun(@mean, Y)];
            
            
            %% ============= plot distributions of distance score
            if compareNtoFRdist ==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname ',' num2str(i) '-' num2str(ii)]);
                xlabel('fr diff (all ngram pairs)');
                
                % --- go through both pairs
                Y = {};
                for pnum = 1:length(indpaircomp)
                    
                    indtmp = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii & ...
                        OUTSTRUCT.All_diffsyl_PairType == indpaircomp(pnum);
                    
                    assert(sum(indtmp) == sum(All_Birdnum==i & All_Neurnum==ii & ...
                        All_PairType==indpaircomp(pnum)), 'old and new dat dont match?');
                    
                    % --------------- extract old data on fr diff
                    frdiff = OUTSTRUCT.All_AbsFRdiff(indtmp);
                    frdiff_NEG = OUTSTRUCT.All_AbsFRdiff_NEG(indtmp);
                    
                    Y{2*pnum-1} = frdiff;
                    Y{2*pnum} = frdiff_NEG;
                    
                    %                 lt_plot_histogram(frdiff, '', 1, 1, '', 1, pcols_pairs{pnum});
                    %                 lt_plot_histogram(frdiff_NEG, '', 1, 1, '', 1, pcols_pairs{pnum});
                end
                
                % === plot
                assert(length(indpaircomp)==2, 'this code assumes comparing 2');
                x = 1:4;
                Xl = {OUTSTRUCT.PairTypesInOrder{indpaircomp(1)}, 'neg', ...
                    OUTSTRUCT.PairTypesInOrder{indpaircomp(2)}, 'neg'};
                lt_plot_MultDist(Y, x, 0, 'k', 0);
                set(gca, 'XTick', x);
                set(gca, 'XTickLabel', Xl);
                rotateXLabels(gca, 45);
            end
        end
    end
    
    
    %% ============================== PLOT SUMMARY ACROSS ALL NEURONS/BIRDS
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for i=1:maxbirds
        birdname = SummaryStruct.birds(i).birdname;
        
        inds = NComp_birdnum==i;
        
        if ~any(inds)
            continue
        end
        
        n = NComp_N(inds,:);
        p = NComp_pval(inds);
        
        % ========== pliot
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) ' mean N, allpairs']);
        ylabel(OUTSTRUCT.PairTypesInOrder{indpaircomp(2)});
        xlabel(OUTSTRUCT.PairTypesInOrder{indpaircomp(1)});
        
        lt_plot_45degScatter(n(:,1)+rand, n(:,2)+rand, 'k', 1);
        
        %             xlim([0.5 max(NComp_N(:))+1)]);
        %             ylim([0.5 max(NComp_N(:))+1)]);
        % --- plot significnat ones
        indtmp = p<0.05;
        plot(n(indtmp,1)+rand, n(indtmp,2)+rand, 'ob');
        
        % ----------- divide by brain region
        nlist = NComp_neurnum(inds);
        bregionlist = {SummaryStruct.birds(i).neurons(nlist).NOTE_Location};
        % --- RA
        indtmp = strcmp(bregionlist, 'RA');
        plot(n(indtmp,1)+rand, n(indtmp,2)+rand, 'xr')
        
    end
    
%     
%     % ==================== [USING LOG10]
%     figcount=1;
%     subplotrows=4;
%     subplotcols=3;
%     fignums_alreadyused=[];
%     hfigs=[];
%     hsplots = [];
%     
%     for i=1:maxbirds
%         birdname = SummaryStruct.birds(i).birdname;
%         
%         inds = NComp_birdnum==i;
%         
%         if ~any(inds)
%             continue
%         end
%         
%         n = log10(NComp_N(inds,:));
%         p = NComp_pval(inds);
%         
%         % ========== pliot
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         title([birdname ',' num2str(i) ' log(mean N), allpairs']);
%         ylabel(OUTSTRUCT.PairTypesInOrder{indpaircomp(2)});
%         xlabel(OUTSTRUCT.PairTypesInOrder{indpaircomp(1)});
%         
%         lt_plot_45degScatter(n(:,1), n(:,2), 'k', 1);
%         
%         xlim([0.5 log10(max(NComp_N(:))+1)]);
%         ylim([0.5 log10(max(NComp_N(:))+1)]);
%         % --- plot significnat ones
%         indtmp = p<0.05;
%         plot(n(indtmp,1), n(indtmp,2), 'or');
%         
%     end

%% ================================== COMPARE LMAN VS. RA (LINE PLOT)

    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for i=1:maxbirds
        birdname = SummaryStruct.birds(i).birdname;
        
        inds = NComp_birdnum==i;
        
        if ~any(inds)
            continue
        end
        
        n = NComp_N(inds,:);
%         p = NComp_pval(inds);
        % ----------- divide by brain region
        nlist = NComp_neurnum(inds);
        bregionlist = {SummaryStruct.birds(i).neurons(nlist).NOTE_Location};
        
        % ========== pliot
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) ' mean N, allpairs']);
        ylabel(OUTSTRUCT.PairTypesInOrder{indpaircomp(2)});
        xlabel('LMAN -- RA');
        ylabel('N (pairtype1 vs pairtype2)');

        % --- LMAN
        x = [1 2];
        indtmp = strcmp(bregionlist, 'LMAN');
        if any(indtmp)
        plot(x, [n(indtmp,1)+rand n(indtmp,2)+rand]', '-b');
        end
        
        % --- RA
        x = [3 4];
        indtmp = strcmp(bregionlist, 'RA');
        if any(indtmp)
        plot(x, [n(indtmp,1)+rand n(indtmp,2)+rand]', '-r');
        end
        
        xlim([0 5]);
        ylim([0 max(n(:))+1]);
    end

    
    %% ========================= COMPARE LMAN VS. RA (DIFFERENCE BTW PAIRS)

    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for i=1:maxbirds
        birdname = SummaryStruct.birds(i).birdname;
        
        inds = NComp_birdnum==i;
        
        if ~any(inds)
            continue
        end
        
        n = NComp_N(inds,:);
%         p = NComp_pval(inds);
        % ----------- divide by brain region
        nlist = NComp_neurnum(inds);
        bregionlist = {SummaryStruct.birds(i).neurons(nlist).NOTE_Location};
        
        % ========== pliot
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) ' mean N, allpairs']);
        xlabel('LMAN -- RA');
        ylabel('N (pairtype2 DIVIDED BY pairtype2)');

        Y = {};
        % --- LMAN
        indtmp = strcmp(bregionlist, 'LMAN');
        Y{1} = n(indtmp,2)./n(indtmp,1);
        % --- RA
        indtmp = strcmp(bregionlist, 'RA');
        Y{2} = n(indtmp,2)./n(indtmp,1);
        
        % ======= plot
        lt_plot_MultDist(Y, [1 2], 1);
        
        xlim([0 3]);
%         ylim([0 max(n(:))+1]);
ylim([0 2]);
line(xlim, [1 1]);
    end

end

%% distribution of sample sizes across neurons for a given bird

indtmp = OUTSTRUCT.All_birdnum==5;
neurnum = OUTSTRUCT.All_neurnum(indtmp);
tabulate(neurnum)


















