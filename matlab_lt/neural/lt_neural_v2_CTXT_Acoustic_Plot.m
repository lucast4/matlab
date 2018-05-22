function lt_neural_v2_CTXT_Acoustic_Plot(analyfname, bregiontoplot)
%% lt 5/19/18 - how does FF variability relate to context coding?
% bregiontoplot = 'LMAN';

%% params

plottext = 0; % neuron number
apos =1;

%% load branch
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

% load([savedir '/CLASSESv2_' analyfname '.mat']);
load([savedir '/ALLBRANCHv2_' analyfname '.mat']);
load([savedir '/SUMMARYv2_' analyfname '.mat']);
try
    load([savedir '/PARAMSv2_' analyfname '.mat']);
catch err
end

%% ======= Filter data (e.g. remove noise, poor labels, etc)

Params.LocationsToKeep = {};
Params.birdstoexclude = {};
Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
Params.durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
Params.durThreshOmega.gappre= [];
Params.durThreshOmega.gappost= [];
Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
% class within branch) is longer than this (sec)
Params.RemoveHandCoded =1 ; % see below

ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params);

%% use first time bin only (can modify)

tt = 1;

%% ========= COMPUILE BY ACTUAL BRANCH IDS (based on regexp str)
DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, ...
    analyfname, tt, 0, [], prms);

disp('MANY EXTRACTION FAILURES BECUASE BRANCH POINTS DOESNT EXIST!! - IS OK')


%% ========= GO THRU ALL BRANCHES, IS THERE RELATIONSHIP BETWEEN DECODE AND DIFF IN MEAN FF?

numbirds = length(DATSTRUCT_BYBRANCH.bird);

AllFFeta2 = [];
AllFFpval = [];
AllDecode = [];
AllDecode_z = [];
AllDecode_p = [];
AllNctxt = [];
AllBird = [];
AllBranch = [];

for i=1:numbirds
    
    numbranch = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb = 1:numbranch
        disp([num2str(i) '-' num2str(bb)]);
        
        DAT = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT;
        regexpthis = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========== get inds for this brain region
        inds = strcmp(DAT.brainregion, bregiontoplot);
        
        % ========== EXTRACT data
        % ---- 1) get difference in FF across contexts for each neuron
        ffstruct = DAT.FFstruct(inds);
        nneur = length(ffstruct);
        EtaSqAll = [];
        PvalAll = [];
        for nn=1:nneur
            %%
            FFeachclass = {};
            FFeachclass_minusmean = {};
            nclass = length(ffstruct(nn).classnum);
            for cc = 1:nclass
                ff = ffstruct(nn).classnum(cc).t_ff(:,2)';
                FFeachclass{cc} = ff;
                FFeachclass_minusmean{cc} = ff - mean(ff);
            end
            
            % == GET AMOUNT OF VARIANCE ACCOUNTED FOR BY CONTEXT
            if (0)
                % ---- total variance
                ffall = cell2mat(FFeachclass);
                vartot = var(ffall);
                
                % ---- variance after subtracting mean
                ffall = cell2mat(FFeachclass_minusmean);
                var_resid = var(ffall);
                
                R2 = (vartot - var_resid)/vartot;
                
                if isnan(R2)
                    % then make sure if becuase FF doesn't exist
                    assert(all(isnan(ffall)), 'why is R2 nan?');
                end
                
                % ========= odl versio, uses anova
            else
                ffall = [];
                classall =[];
                nclass = length(ffstruct(nn).classnum);
                for cc = 1:nclass
                    ffthis = ffstruct(nn).classnum(cc).t_ff(:,2)';
                    ffall = [ffall ffthis];
                    classall = [classall cc*ones(size(ffthis))];
                end
                
                [p, anovatbl] = anova1(ffall, classall, 'off');
                
                ss_effect = anovatbl{2,2};
                ss_total = anovatbl{4,2};
                df_effect = anovatbl{2,3};
                ms_error = anovatbl{3,4};
                
                
                % calculate omega squared
                numerator = ss_effect - df_effect*ms_error;
                denominator = ss_total + ms_error;
                omega2 = numerator/denominator;
                
                % calculate eta squared
                R2 = ss_effect / ss_total;
            end
                EtaSqAll = [EtaSqAll; R2];
                PvalAll = [PvalAll; p];
            
        end
        
        % ---------- FR DAT
        frdat = DAT.frdat(inds);
        
        % --------  2) Get decode
        ydecode_prem = [];
        ydecode_prem_z = [];
        for k = find(inds)'
            sts = lt_neural_ConfMatStats(DAT.PREMOTORDECODE_struct(k).ConfMatAll_DAT);
            ydecode_prem = [ydecode_prem; mean([sts.F1])];
            
            % ---- get zscore rel shuffle
            stsneg = lt_neural_ConfMatStats(DAT.PREMOTORDECODE_struct(k).ConfMatAll_NEG(1:200));
            
            ydecode_prem_z = [ydecode_prem_z; ...
                (mean([sts.F1]) - mean([stsneg.F1]))/std([stsneg.F1])];
            
        end
        ydecode_pval = DAT.PREMOTORDECODE_pval(inds);
        
        % ---------- GET DECODE (zscore rel base)
        
        
        % --------- 3) Number of contexts
        Nctxt = cellfun(@length, {DAT.frdat(inds).classnum})';
        
        
        % ============================== OUTPUT
        AllFFeta2 = [AllFFeta2; EtaSqAll];
        AllFFpval = [AllFFpval; PvalAll];
        AllDecode = [AllDecode; ydecode_prem];
        AllDecode_z = [AllDecode_z; ydecode_prem_z];
        AllNctxt = [AllNctxt; Nctxt];
        AllBird = [AllBird; i*ones(size(Nctxt))];
        AllBranch = [AllBranch; bb*ones(size(Nctxt))];
        AllDecode_p = [AllDecode_p; ydecode_pval];
        
    end
end



%% ================= [PLOT] decode vs. FF R2, as function of nctxt

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

%
% for i=1:maxNct
%
%     for bb=1:maxBirds
%
%     inds = AllNctxt==i & AllBird==bb;
%
%     if ~any(~isnan(AllFFeta2(inds)))
%         continue
%     end
%
%     % ===================== COLLECT
%     ffeta = AllFFeta2(inds);
%     ydecode = AllDecode(inds);
%
%     % ================== PLOT
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     birdname = SummaryStruct.birds(bb).birdname;
%
%     title(['nctxt: ' num2str(i) '[' birdname  ']']);
%     xlabel('ff diff (anova r2)');
%     ylabel('neural decode');
%     plot(ffeta, ydecode, 'oR');
% %     xlim([-0.1 1]);
% %     ylim([-0.1 1]);
% %     lt_plot_zeroline;
% %     lt_plot_zeroline_vert;
%     end
% end
%


pcolors = lt_make_plot_colors(maxBirds, 0,0);

for bb=1:maxBirds
    pcol = pcolors{bb};
    
    for i=1:maxNct
        
        pcol = [rand rand rand];
        
        for ii = 1:maxBranchnum
            
            inds = AllNctxt==i & AllBird==bb & AllBranch==ii;
            
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ffeta = AllFFeta2(inds);
            ydecode = AllDecode(inds);
            
            % ================== PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            birdname = SummaryStruct.birds(bb).birdname;
            branchname = DATSTRUCT_BYBRANCH.bird(bb).branchID(ii).regexpstr{1};
            
            title(['nctxt: ' num2str(i) '[' birdname  ']' branchname]);
            xlabel('ff diff (anova r2)');
            ylabel('neural decode');
            plot(ffeta, ydecode, 'o', 'Color',pcol);
            xlim([-0.1 1]);
            ylim([-0.1 1]);
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
        end
    end
end



%% =============== [PLOT] DECODE VS. FF R2, in combined plot


maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



pcolors = lt_make_plot_colors(maxBirds, 0,0);

for bb=1:maxBirds
    
    for i=1:maxNct
        
        inds = AllNctxt==i & AllBird==bb;
        
        if ~any(~isnan(AllFFeta2(inds)))
            continue
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        birdname = SummaryStruct.birds(bb).birdname;
        
        title(['nctxt: ' num2str(i) '[' birdname  ']']);
        xlabel('ff diff (anova r2)');
        ylabel('neural decode');
        
        for ii = 1:maxBranchnum
            
            inds = AllNctxt==i & AllBird==bb & AllBranch==ii;
            
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ffeta = AllFFeta2(inds);
            ydecode = AllDecode(inds);
            
            % ================== PLOT
            pcol = [rand rand rand];
            plot(ffeta, ydecode, 'o', 'Color',pcol);
            xlim([-0.1 1]);
            ylim([-0.1 1]);
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
        end
    end
end


%% =============== [PLOT] DECODE VS. FF R2, in combined plot
% ZSCORE

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



pcolors = lt_make_plot_colors(maxBirds, 0,0);

for bb=1:maxBirds
    
    for i=1:maxNct
        
        inds = AllNctxt==i & AllBird==bb;
        
        if ~any(~isnan(AllFFeta2(inds)))
            continue
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        birdname = SummaryStruct.birds(bb).birdname;
        
        title(['nctxt: ' num2str(i) '[' birdname  ']']);
        xlabel('ff diff (anova r2)');
        ylabel('neural decode (zscore)');
        
        for ii = 1:maxBranchnum
            
            inds = AllNctxt==i & AllBird==bb & AllBranch==ii;
            
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ffeta = AllFFeta2(inds);
            ydecode = AllDecode_z(inds);
            ydecode_p = AllDecode_p(inds);
            
            % ================== PLOT
            pcol = [rand rand rand];
            
            % --- go thrue each and plot solid if decode is sig
            for j=1:length(ydecode)
                if ydecode_p(j)<0.05
            lt_plot(ffeta(j), ydecode(j), {'Color',pcol});        
                else
            plot(ffeta(j), ydecode(j), 'o', 'Color',pcol);        
                end
            end
            xlim([-0.1 1]);
            ylim([-0.5 max(AllDecode_z)]);
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
        end
    end
end



%% =============== [PLOT] DECODE VS. FF R2, in combined plot
% TAKE MEAN ACROSS ALL NEURONS FOR A GIVEN BRANCH

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

for i=1:maxNct
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    birdname = SummaryStruct.birds(bb).birdname;
    
    title(['nctxt: ' num2str(i) '[allbirds]']);
    xlabel('ff diff (anova r2)');
    ylabel('neural decode');
    
    for bb=1:maxBirds
        
        %         inds = AllNctxt==i & AllBird==bb;
        %
        %         if ~any(~isnan(AllFFeta2(inds)))
        %             continue
        %         end
        %
        pcol = pcolsbirds{bb};
        
        for ii = 1:maxBranchnum
            
            inds = AllNctxt==i & AllBird==bb & AllBranch==ii;
            
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ffeta = AllFFeta2(inds);
            ydecode = AllDecode(inds);
           
            
            % --------- take means
            ffeta_sem = lt_sem(ffeta);
            ffeta = mean(ffeta);
            ydecode_sem = lt_sem(ydecode);
            ydecode = mean(ydecode);
            
            % ================== PLOT
            lt_plot(ffeta, ydecode, {'Color',pcol});
            
            line(ffeta+[-ffeta_sem ffeta_sem], [ydecode ydecode], 'Color', pcol);
            line([ffeta ffeta], ydecode+[-ydecode_sem ydecode_sem], 'Color', pcol);
            
            xlim([-0.1 1]);
            ylim([-0.1 1]);
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
        end
    end
end





%% =============== [PLOT] DECODE VS. FF R2, in combined plot
% TAKE MEAN ACROSS ALL NEURONS FOR A GIVEN BRANCH 

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

for i=1:maxNct
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    birdname = SummaryStruct.birds(bb).birdname;
    
    title(['nctxt: ' num2str(i) '[allbirds]']);
    xlabel('ff diff (anova r2)');
    ylabel('neural decode (zscore)');
    
    for bb=1:maxBirds
        
        %         inds = AllNctxt==i & AllBird==bb;
        %
        %         if ~any(~isnan(AllFFeta2(inds)))
        %             continue
        %         end
        %
        pcol = pcolsbirds{bb};
        
        for ii = 1:maxBranchnum
            
            inds = AllNctxt==i & AllBird==bb & AllBranch==ii;
            
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ffeta = AllFFeta2(inds);
            ydecode = AllDecode_z(inds);
            ydecode_p = AllDecode_p(inds) <0.05;
            
            % --------- take means
            ffeta_sem = lt_sem(ffeta);
            ffeta = mean(ffeta);
            ydecode_sem = lt_sem(ydecode);
            ydecode = mean(ydecode);
            
            if sum(ydecode_p)>=max([1 floor(length(ydecode_p)/2)])
                % then at least half significant
                plotsig = 1;
            else
                plotsig = 0;
            end
                            
            % ================== PLOT
            if plotsig==1
            lt_plot(ffeta, ydecode, {'Color',pcol});
            else
                plot(ffeta, ydecode, 'o', 'Color',pcol);
            end
            line(ffeta+[-ffeta_sem ffeta_sem], [ydecode ydecode], 'Color', pcol);
            line([ffeta ffeta], ydecode+[-ydecode_sem ydecode_sem], 'Color', pcol);
            
            xlim([-0.1 1]);
            ylim([-0.5 max(AllDecode_z)+0.1]);
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
        end
    end
end





%% =================== [PLOT, SEPARATE BY IF FF DIFFER SIGNIFICANTLY]
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);

pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

    
    for bb=1:maxBirds
        
        inds = AllBird==bb;
        
        if ~any(~isnan(AllFFeta2(inds)))
            continue
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        birdname = SummaryStruct.birds(bb).birdname;
        
        title(['[' birdname ']']);
        xlabel('numcases (open = sign, ff diff)');
        ylabel('numcases (significant decode)');
        pcol = pcolsbirds{bb};
        
        for ii = 1:maxBranchnum
            
            inds = AllBird==bb & AllBranch==ii;
            pcol = [rand rand rand];
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ff_p = AllFFpval(inds);
            ydecode_p = AllDecode_p(inds);
            
            
            % ==================== separate by whether ff is significantly
            % diff
            
            
            % ------------ sig diff ff
            indtmp = ff_p>=0.05;
            
            numcases = sum(indtmp);
            numcases_sigdecode = sum(ydecode_p(indtmp)<0.05);
            x1 = numcases+0.4*rand-0.2;
            y1 = numcases_sigdecode+0.4*rand-0.2;
            lt_plot(x1, y1, {'Color', pcol});
            
            % ------------ sig diff ff
            indtmp = ff_p<0.05;
            
            numcases = sum(indtmp);
            numcases_sigdecode = sum(ydecode_p(indtmp)<0.05);
            x2 = numcases+0.4*rand-0.2;
            y2 = numcases_sigdecode+0.4*rand-0.2;
            
            plot(x2, y2, 'o', 'Color', pcol);
            
            % --- connect by line
            line([x1 x2], [y1 y2], 'Color', pcol);
            lt_plot_makesquare_plot45line(gca, 'k', -0.5);
        end
    end

%% =================== [PLOT, SEPARATE BY IF FF DIFFER SIGNIFICANTLY]
% ONE PLOT, *eacg bird diff color
lt_figure; hold on;
        title(['[allbird]']);
        xlabel('numcases (open = sign, ff diff)');
        ylabel('numcases (significant decode)');

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);

pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

    
    for bb=1:maxBirds
        
        inds = AllBird==bb;
        
        if ~any(~isnan(AllFFeta2(inds)))
            continue
        end
        
        pcol = pcolsbirds{bb};
        
            
            inds = AllBird==bb 
            if ~any(~isnan(AllFFeta2(inds)))
                continue
            end
            
            % ===================== COLLECT
            ff_p = AllFFpval(inds);
            ydecode_p = AllDecode_p(inds);
            
            
            % ==================== separate by whether ff is significantly
            % diff
            
            
            % ------------ not sig diff ff
            indtmp = ff_p>=0.05;
            
            numcases = sum(indtmp);
            numcases_sigdecode = sum(ydecode_p(indtmp)<0.05);
            x1 = numcases+0.4*rand-0.2;
            y1 = numcases_sigdecode+0.4*rand-0.2;
            lt_plot(x1, y1, {'Color', pcol});
            
            % ------------ sig diff ff
            indtmp = ff_p<0.05;
            
            numcases = sum(indtmp);
            numcases_sigdecode = sum(ydecode_p(indtmp)<0.05);
            x2 = numcases+0.4*rand-0.2;
            y2 = numcases_sigdecode+0.4*rand-0.2;
            
            plot(x2, y2, 'o', 'Color', pcol);
            
            % --- connect by line
            line([x1 x2], [y1 y2], 'Color', pcol);
            lt_plot_makesquare_plot45line(gca, 'k', -0.5);
        end


%% =================== [PLOT, DECODE PVAL AS FUNCTION OF FF PVAL]
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);

pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

for bb=1:maxBirds
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    birdname = SummaryStruct.birds(bb).birdname;
    title(['[' birdname ']']);
    xlabel('FF (anova pval)');
    ylabel('Decode (pval)');
        
    for ii = 1:maxBranchnum
        
        inds = AllBird==bb & AllBranch==ii;
        if ~any(~isnan(AllFFeta2(inds)))
            continue
        end
        
        pcol = 0.7*[rand rand rand];
        
        % ===================== COLLECT
        ff_p = log10(AllFFpval(inds));
        % -- lower lim for p val at 0.001
        ff_p(ff_p<-3) = -3;
        ydecode_p = log10(AllDecode_p(inds));
        
        plot(ff_p, ydecode_p, 'o', 'Color', pcol);
                
        % -----
        line(xlim, log10([0.05 0.05]));
        line(log10([0.05 0.05]), ylim);
%         lt_plot_makesquare_plot45line(gca, 'k', -0.5);
        xlim([-3.2 0]);
        ylim([-3.2 0]);
    end
end




%% =================== [PLOT, DECODE ZSCORE AS FUNCTION OF FF R2]
% GOOD: both effect size and significance ...

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

maxBirds = max(AllBird);
maxNct = max(AllNctxt);
maxBranchnum = max(AllBranch);

pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

for bb=1:maxBirds
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    birdname = SummaryStruct.birds(bb).birdname;
    title(['[' birdname ']']);
    xlabel('FF (anova R2)[sig=circle]');
    ylabel('Decode (z-score) [solid=sig]');
        
    for ii = 1:maxBranchnum
        
        inds = AllBird==bb & AllBranch==ii;
        if ~any(~isnan(AllFFeta2(inds)))
            continue
        end
        
        pcol = 0.7*[rand rand rand];
        
        % ===================== COLLECT
        ff_r2 = AllFFeta2(inds);
        ff_p = AllFFpval(inds);
        ydecode_z = AllDecode_z(inds);
        ydecode_p = AllDecode_p(inds);
       
        
        % ---- significant decodes/not sig FF
        indtmp = ydecode_p<0.05 & ff_p>=0.05;
        lt_plot(ff_r2(indtmp), ydecode_z(indtmp), {'Color', pcol, 'Marker', 'p'});
                
        % --- signfiicant decode/significant FF
        indtmp = ydecode_p<0.05 & ff_p<0.05;
        lt_plot(ff_r2(indtmp), ydecode_z(indtmp), {'Color', pcol, 'Marker', 'o', 'MarkerSize', 5});
        
        % ---- notsignificant decodes/not sig FF
        indtmp = ydecode_p>=0.05 & ff_p>=0.05;
        plot(ff_r2(indtmp), ydecode_z(indtmp), 'p', 'Color', pcol);
                
        % --- not signfiicant decode/significant FF
        indtmp = ydecode_p>=0.05 & ff_p<0.05;
        plot(ff_r2(indtmp), ydecode_z(indtmp), 'o', 'Color', pcol);

        % -----
        xlim([-0.1 1.1]);
        ylim([min(AllDecode_z)-0.1 max(AllDecode_z)+0.1]);
        lt_plot_zeroline;
        
%         line(xlim, log10([0.05 0.05]));
%         line(log10([0.05 0.05]), ylim);
% %         lt_plot_makesquare_plot45line(gca, 'k', -0.5);
%         xlim([-3.2 0]);
%         ylim([-3.2 0]);
    end
end





%% =================== [PLOT, DECODE ZSCORE AS FUNCTION OF FF R2]
% GOOD: both effect size and significance ... [all combined]
lt_figure; hold on;
title(['[all birds]']);
xlabel('FF (anova R2)[sig=circle]');
ylabel('Decode (z-score) [solid=sig]');

maxBirds = max(AllBird);
maxBranchnum = max(AllBranch);
pcolsbirds = lt_make_plot_colors(maxBirds, 0,0);

for bb=1:maxBirds
    
    inds = AllBird==bb;
    if ~any(~isnan(AllFFeta2(inds)))
        continue
    end
    
    pcol = pcolsbirds{bb};
    
    % ===================== COLLECT
    ff_r2 = AllFFeta2(inds);
    ff_p = AllFFpval(inds);
    ydecode_z = AllDecode_z(inds);
    ydecode_p = AllDecode_p(inds);
    
    
    % ---- significant decodes/not sig FF
    indtmp = ydecode_p<0.05 & ff_p>=0.05;
    lt_plot(ff_r2(indtmp), ydecode_z(indtmp), {'Color', pcol, 'Marker', 'p'});
    
    % --- signfiicant decode/significant FF
    indtmp = ydecode_p<0.05 & ff_p<0.05;
    lt_plot(ff_r2(indtmp), ydecode_z(indtmp), {'Color', pcol, 'Marker', 'o', 'MarkerSize', 5});
    
    % ---- notsignificant decodes/not sig FF
    indtmp = ydecode_p>=0.05 & ff_p>=0.05;
    plot(ff_r2(indtmp), ydecode_z(indtmp), 'p', 'Color', pcol);
    
    % --- not signfiicant decode/significant FF
    indtmp = ydecode_p>=0.05 & ff_p<0.05;
    plot(ff_r2(indtmp), ydecode_z(indtmp), 'o', 'Color', pcol);
    
end

    % -----
    xlim([-0.1 1.1]);
    ylim([min(AllDecode_z)-0.1 max(AllDecode_z)+0.1]);
    lt_plot_zeroline;
    


%% ======  [PROB OF SIGNIFICANT DECODE, DEPEND ON FF DIFF SIGNIFANT?];

lt_figure; hold on;
title(['[all birds]']);
xlabel('FFdiff [not sig -- sig]');
ylabel('frac cases sig decode');


% ===================== COLLECT
inds = ~isnan(AllFFeta2); % --- only keep cases with FF
% ff_r2 = AllFFeta2(inds);
ff_p = AllFFpval(inds);
% ydecode_z = AllDecode_z(inds);
ydecode_p = AllDecode_p(inds);

% ============== collect fractions
Y = [];
% ---------------- not signinfinct FF
indtmp = ff_p>=0.05;
nsig_decode = sum(ydecode_p(indtmp)<0.05);
ntot_decode = sum(indtmp);

Y = [Y nsig_decode/ntot_decode];

% ---------------- signinfinct FF
indtmp = ff_p<0.05;
nsig_decode = sum(ydecode_p(indtmp)<0.05);
ntot_decode = sum(indtmp);

Y = [Y nsig_decode/ntot_decode];


% -------- plot
X = [1 2];

plot(X, Y, '-ok');
xlim([0 3]);
lt_plot_zeroline;






