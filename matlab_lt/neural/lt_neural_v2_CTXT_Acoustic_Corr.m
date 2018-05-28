function lt_neural_v2_CTXT_Acoustic_Corr(analyfname, bregiontoplot, doshuffmany)

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


%% ========= determine premotor windows

frinds = 1:1000*(prms.motifpredur + prms.motifpostdur)-1;

%%

AllClassnum = {};
AllPitch = {};
AllFRbinned = {};
AllFRmeans = {};
AllBirdnum = [];
AllBranchnum = [];
AllNeurnum = [];
AllNeurnum_fakeInd = [];

numbirds = length(DATSTRUCT_BYBRANCH.bird);
for i =1:numbirds
    numbranch = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    for bb=1:numbranch
        
        DAT = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT;
        indsneur = find(strcmp(DAT.brainregion, bregiontoplot));
        
        
        for nn=indsneur'
            
            nclass = length(DAT.frdat(nn).classnum);
            
            
            % ================ pitch and class and frmat
            pitchall = [];
            classall = [];
            frmatall = [];
            for cc=1:nclass
                ff = DAT.FFstruct(nn).classnum(cc).t_ff(:,2)';
                frmat = DAT.frdat(nn).classnum(cc).frmat(frinds,:); % take frinds, so same timebins.
                
                pitchall = [pitchall ff];
                classall = [classall cc*ones(size(ff))];
                frmatall = [frmatall frmat];
            end
            
            % ------ take premotor window for frmat (for frmat, determine which time inds to take (i.e.
            % premotor window)
            premwind = DAT.PREMOTORDECODE_struct(nn).window_relonset;
            xwind = prms.motifpredur + premwind; % window, in sec
            
            % ---- if window too short, lengthen to 10ms.
            if xwind(2)-xwind(1)<0.01
                xwind(2) = xwind(1)+0.01;
            end
            
            xtimes = lt_neural_QUICK_XfromFRmat(frmatall); % time of bins
            
            indsgood = xtimes>=xwind(1) & xtimes<xwind(2);
            frmatall = frmatall(indsgood, :);
            xtimes = xtimes(indsgood);
            
            % ------ convert frmatall to fr bins
            TrimDown = 1;
            frbinsize = 0.008;
            [frbinned, xtimes] = lt_neural_v2_QUICK_binFR(frmatall', xtimes, frbinsize, TrimDown);
            
            % ----- get one value for mean fr in window (i.e. spike count)
            frmeans = mean(frmatall,1);
            
            
            % ======================= OUTPUT
            AllClassnum = [AllClassnum; classall];
            AllPitch = [AllPitch; pitchall];
            AllFRbinned = [AllFRbinned; frbinned];
            AllFRmeans = [AllFRmeans; frmeans];
            AllBirdnum = [AllBirdnum; i];
            AllBranchnum = [AllBranchnum; bb];
            
            AllNeurnum = [AllNeurnum; DAT.IndsOriginal_neuron(nn)];
            AllNeurnum_fakeInd = [AllNeurnum_fakeInd; nn];
            
        end
    end
end

%%
maxbranch = max(AllBranchnum);
maxneur = max(AllNeurnum);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [CORR - POOL ACROSS CONTEXTS]
%% ========== is there significant correlation between Frate and pitch?
% [SINGLE INSTANCE SHUFFLES]

poolOverCtxt = 1;

doshuff = 0;
[rhoall, rhoCIall, rhoPall, Nsampall] = fn_getallRho(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff, poolOverCtxt);

% ========= do one shuffle cycle [within each class]
doshuff = 1;
[rhoallSHUFF, rhoCIallSHUFF, rhoPallSHUFF, NsampallSHUFF] = fn_getallRho(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff, poolOverCtxt);


% ========================== PLOT
lt_figure; hold on;

% ---------------- 1) ALL RHO
lt_subplot(3,2,1); hold on;
title('all rho');
xlabel('rho');
[~, xcenters] = lt_plot_histogram(rhoall, '', 1, 0, '', 1, 'k');
% ---- overlay shuffle
lt_plot_histogram(rhoallSHUFF, xcenters, 1, 0, '', 1, 'r');
% ---- more pos or negative corrs on average?
p = signrank(rhoall, rhoallSHUFF);
lt_plot_pvalue(p, 'srank', 1);

% ---------------- 2) absolute value of corr (dat vs. shuff);
lt_subplot(3,2,3); hold on;
xlabel('abs(rho) (dat)');
ylabel('abs(rho) (SHUFF)');
x = abs(rhoall);
y = abs(rhoallSHUFF);
plot(x, y, 'o');
xlim([-0.1 1.1]);
ylim([-0.1 1.1]);
lt_plot_makesquare_plot45line(gca, 'b', -0.1);
p = signrank(x, y);
lt_plot_pvalue(p, 'srank', 1);

% ---------------- abs value, but paired lines
lt_subplot(3,2,5); hold on;
xlabel('dat -- shuff');
ylabel('abs(rho)');
x = [1 2];
y = [abs(rhoall) abs(rhoallSHUFF)];
plot(x, y', '-o', 'Color', [0.7 0.7 0.7]);
lt_plot(x+0.2, mean(y), {'Errors', lt_sem(y)});
xlim([0 3]);
lt_plot_zeroline;
p = signrank(y(:,1), y(:,2));
lt_plot_pvalue(p, 'srank', 1);



% -------------- 2) RHO (WITH CI) grouped by sample size)
lt_subplot(3,2,2); hold on;
xlabel('rho(CI)');
ylabel('log2(N (samp))');
% --- significant corr
indtmp = rhoPall<0.05;
lt_plot(rhoall(indtmp), log2(Nsampall(indtmp)), {'Color', 'k'});
% --- not sig
indtmp = rhoPall>=0.05;
plot(rhoall(indtmp), log2(Nsampall(indtmp)), 'o', 'Color', 'k');

% ---------------- same, but for shuffle
lt_subplot(3,2,4); hold on;
title('SHUFF', 'Color', 'r');
xlabel('rho(CI)');
ylabel('log2(N (samp))');
% --- significant corr
indtmp = rhoPallSHUFF<0.05;
lt_plot(rhoallSHUFF(indtmp), log2(NsampallSHUFF(indtmp)), {'Color', 'r'});
% --- not sig
indtmp = rhoPallSHUFF>=0.05;
plot(rhoallSHUFF(indtmp), log2(NsampallSHUFF(indtmp)), 'o', 'Color', 'r');


% ============ percent of cases significatn (extrapolate based on sampsize)
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
xlabel('sample size');
ylabel('frac cases sig');
for n=unique(Nsampall)'
    
    indtmp = Nsampall==n;
    nsig = sum(rhoPall(indtmp)<0.05)/length(rhoPall(indtmp));
    plot(n, nsig, 'ok');
end


lt_subplot(3,2,2); hold on;
title('divide data based on samp size');
xlabel('median sampel size');
ylabel('frac sig');
sampedges = prctile(Nsampall, [25 50 75 100]);
sampedges = [0 sampedges];
for i=1:length(sampedges)-1
    
    indtmp = Nsampall>sampedges(i) & Nsampall<=sampedges(i+1);
    
    nsig = sum(rhoPall(indtmp)<0.05)/length(rhoPall(indtmp));
    sampmedian = median(Nsampall(indtmp));
    
    plot(sampmedian, nsig, '-ok');
end

lt_plot_zeroline;
lt_plot_zeroline_vert;

%% ============== [MULTIPLE SHUFFLES]
if doshuffmany == 1
    % ================== COLLLECT DATA
    doshuff = 0;
    poolOverCtxt=1;
    [rhoallDAT, rhoCIallDAT, rhoPallDAT, NsampallDAT] = fn_getallRho(AllBirdnum, AllBranchnum, ...
        AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff, poolOverCtxt);
    
    % ================= COLLECT MULTIPLE SHUFFLES
    Nshuffs = 1000;
    AbsRhoMedian_ShuffMult = nan(Nshuffs, 1);
    RhoMedian_ShuffMult = nan(Nshuffs, 1);
    Nsign_ShuffMult = nan(Nshuffs,1);
    NsignRight_ShuffMult = nan(Nshuffs,1);
    NsignLeft_ShuffMult = nan(Nshuffs,1);
    for i=1:Nshuffs
        disp(num2str(i));
        doshuff = 1;
        [rhoallSHUFF, ~, rhoPallSHUFF] = fn_getallRho(AllBirdnum, AllBranchnum, ...
            AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff, poolOverCtxt);
        
        % ==== collect median abs value of rho
        RhoMedian_ShuffMult(i) = median(rhoallSHUFF);
        AbsRhoMedian_ShuffMult(i) = median(abs(rhoallSHUFF));
        Nsign_ShuffMult(i) = sum(rhoPallSHUFF<0.05);
        NsignRight_ShuffMult(i) = sum(rhoPallSHUFF<0.05 & rhoallSHUFF>0);
        NsignLeft_ShuffMult(i) = sum(rhoPallSHUFF<0.05 & rhoallSHUFF<0);
    end
    
    %% ========================= [PLOTS]
    lt_figure; hold on;
    
    % ------------ 1) median of abs corr
    lt_subplot(3,2,1); hold on;
    title('median of abs corr');
    
    lt_plot_histogram(AbsRhoMedian_ShuffMult);
    ydat = median(abs(rhoallDAT));
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(AbsRhoMedian_ShuffMult>=ydat))./length(AbsRhoMedian_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    
    
    % --------------- median of corr
    lt_subplot(3,2,5); hold on;
    title('median of corr');
    
    lt_plot_histogram(RhoMedian_ShuffMult);
    ydat = median(rhoallDAT);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(abs(RhoMedian_ShuffMult)>=abs(ydat)))./length(RhoMedian_ShuffMult);
    lt_plot_pvalue(p, 'two tailed', 1);
    
    
    % ------------- 2) number significant
    lt_subplot(3,2,2); hold on;
    title('number sig cases');
    xcenters = 0:1:sum(rhoPallDAT<0.05)+1;
    lt_plot_histogram(Nsign_ShuffMult, xcenters);
    ydat = sum(rhoPallDAT<0.05);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(Nsign_ShuffMult>=ydat))/length(Nsign_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    lt_plot_annotation(1, ['N=' num2str(length(rhoallDAT))], 'b');
    
    % ------------ 3) number significant (right);
    lt_subplot(3,2,3); hold on;
    title('number sig cases (positive)');
    
    lt_plot_histogram(NsignRight_ShuffMult, xcenters);
    ydat = sum(rhoPallDAT<0.05 & rhoallDAT>0);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(NsignRight_ShuffMult>=ydat))/length(NsignRight_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    
    
    % ------------ 3) number significant (left);
    lt_subplot(3,2,4); hold on;
    title('number sig cases (negative)');
    
    lt_plot_histogram(NsignLeft_ShuffMult, xcenters);
    ydat = sum(rhoPallDAT<0.05 & rhoallDAT<0);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(NsignLeft_ShuffMult>=ydat))/length(NsignLeft_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [CORR - TREAT EACH CONTEXT AS SYL]
%% ========== is there significant correlation between Frate and pitch?
% ===== for each class, do separate regression. get distribution of
% regression coefficients
% [SINGLE INSTANCE SHUFFLES]

doshuff = 0;
[rhoall, rhoCIall, rhoPall, Nsampall] = fn_getallRho(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff);

% ========= do one shuffle cycle [within each class]
doshuff = 1;
[rhoallSHUFF, rhoCIallSHUFF, rhoPallSHUFF, NsampallSHUFF] = fn_getallRho(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff);


% ========================== PLOT
lt_figure; hold on;

% ---------------- 1) ALL RHO
lt_subplot(3,2,1); hold on;
title('all rho');
xlabel('rho');
[~, xcenters] = lt_plot_histogram(rhoall, '', 1, 0, '', 1, 'k');
% ---- overlay shuffle
lt_plot_histogram(rhoallSHUFF, xcenters, 1, 0, '', 1, 'r');
% ---- more pos or negative corrs on average?
p = signrank(rhoall, rhoallSHUFF);
lt_plot_pvalue(p, 'srank', 1);

% ---------------- 2) absolute value of corr (dat vs. shuff);
lt_subplot(3,2,3); hold on;
xlabel('abs(rho) (dat)');
ylabel('abs(rho) (SHUFF)');
x = abs(rhoall);
y = abs(rhoallSHUFF);
plot(x, y, 'o');
xlim([-0.1 1.1]);
ylim([-0.1 1.1]);
lt_plot_makesquare_plot45line(gca, 'b', -0.1);
p = signrank(x, y);
lt_plot_pvalue(p, 'srank', 1);

% ---------------- abs value, but paired lines
lt_subplot(3,2,5); hold on;
xlabel('dat -- shuff');
ylabel('abs(rho)');
x = [1 2];
y = [abs(rhoall) abs(rhoallSHUFF)];
plot(x, y', '-o', 'Color', [0.7 0.7 0.7]);
lt_plot(x+0.2, mean(y), {'Errors', lt_sem(y)});
xlim([0 3]);
lt_plot_zeroline;
p = signrank(y(:,1), y(:,2));
lt_plot_pvalue(p, 'srank', 1);



% -------------- 2) RHO (WITH CI) grouped by sample size)
lt_subplot(3,2,2); hold on;
xlabel('rho(CI)');
ylabel('log2(N (samp))');
% --- significant corr
indtmp = rhoPall<0.05;
lt_plot(rhoall(indtmp), log2(Nsampall(indtmp)), {'Color', 'k'});
% --- not sig
indtmp = rhoPall>=0.05;
plot(rhoall(indtmp), log2(Nsampall(indtmp)), 'o', 'Color', 'k');

% ---------------- same, but for shuffle
lt_subplot(3,2,4); hold on;
title('SHUFF', 'Color', 'r');
xlabel('rho(CI)');
ylabel('log2(N (samp))');
% --- significant corr
indtmp = rhoPallSHUFF<0.05;
lt_plot(rhoallSHUFF(indtmp), log2(NsampallSHUFF(indtmp)), {'Color', 'r'});
% --- not sig
indtmp = rhoPallSHUFF>=0.05;
plot(rhoallSHUFF(indtmp), log2(NsampallSHUFF(indtmp)), 'o', 'Color', 'r');


% ============ percent of cases significatn (extrapolate based on sampsize)
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
xlabel('sample size');
ylabel('frac cases sig');
for n=unique(Nsampall)'
    
    indtmp = Nsampall==n;
    nsig = sum(rhoPall(indtmp)<0.05)/length(rhoPall(indtmp));
    plot(n, nsig, 'ok');
end


lt_subplot(3,2,2); hold on;
title('divide data based on samp size');
xlabel('median sampel size');
ylabel('frac sig');
sampedges = prctile(Nsampall, [25 50 75 100]);
sampedges = [0 sampedges];
for i=1:length(sampedges)-1
    
    indtmp = Nsampall>sampedges(i) & Nsampall<=sampedges(i+1);
    
    nsig = sum(rhoPall(indtmp)<0.05)/length(rhoPall(indtmp));
    sampmedian = median(Nsampall(indtmp));
    
    plot(sampmedian, nsig, '-ok');
end

lt_plot_zeroline;
lt_plot_zeroline_vert;


%% ============== [MULTIPLE SHUFFLES]
if doshuffmany == 1
    % ================== COLLLECT DATA
    doshuff = 0;
    [rhoallDAT, rhoCIallDAT, rhoPallDAT, NsampallDAT] = fn_getallRho(AllBirdnum, AllBranchnum, ...
        AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff);
    
    % ================= COLLECT MULTIPLE SHUFFLES
    Nshuffs = 1000;
    AbsRhoMedian_ShuffMult = nan(Nshuffs, 1);
    RhoMedian_ShuffMult = nan(Nshuffs, 1);
    Nsign_ShuffMult = nan(Nshuffs,1);
    NsignRight_ShuffMult = nan(Nshuffs,1);
    NsignLeft_ShuffMult = nan(Nshuffs,1);
    for i=1:Nshuffs
        disp(num2str(i));
        doshuff = 1;
        [rhoallSHUFF, ~, rhoPallSHUFF] = fn_getallRho(AllBirdnum, AllBranchnum, ...
            AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff);
        
        % ==== collect median abs value of rho
        RhoMedian_ShuffMult(i) = median(rhoallSHUFF);
        AbsRhoMedian_ShuffMult(i) = median(abs(rhoallSHUFF));
        Nsign_ShuffMult(i) = sum(rhoPallSHUFF<0.05);
        NsignRight_ShuffMult(i) = sum(rhoPallSHUFF<0.05 & rhoallSHUFF>0);
        NsignLeft_ShuffMult(i) = sum(rhoPallSHUFF<0.05 & rhoallSHUFF<0);
    end
    %% ========================= [PLOTS]
    lt_figure; hold on;
    
    % ------------ 1) median of abs corr
    lt_subplot(3,2,1); hold on;
    title('median of abs corr');
    
    lt_plot_histogram(AbsRhoMedian_ShuffMult);
    ydat = median(abs(rhoallDAT));
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(AbsRhoMedian_ShuffMult>=ydat))./length(AbsRhoMedian_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    
    
    % --------------- median of corr
    lt_subplot(3,2,5); hold on;
    title('median of corr');
    
    lt_plot_histogram(RhoMedian_ShuffMult);
    ydat = median(rhoallDAT);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(abs(RhoMedian_ShuffMult)>=abs(ydat)))./length(RhoMedian_ShuffMult);
    lt_plot_pvalue(p, 'two tailed', 1);
    
    
    % ------------- 2) number significant
    lt_subplot(3,2,2); hold on;
    title('number sig cases');
    xcenters = 0:1:sum(rhoPallDAT<0.05)+1;
    lt_plot_histogram(Nsign_ShuffMult, xcenters);
    ydat = sum(rhoPallDAT<0.05);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(Nsign_ShuffMult>=ydat))/length(Nsign_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    lt_plot_annotation(1, ['N=' num2str(length(rhoallDAT))], 'b');
    
    % ------------ 3) number significant (right);
    lt_subplot(3,2,3); hold on;
    title('number sig cases (positive)');
    
    lt_plot_histogram(NsignRight_ShuffMult, xcenters);
    ydat = sum(rhoPallDAT<0.05 & rhoallDAT>0);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(NsignRight_ShuffMult>=ydat))/length(NsignRight_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    
    
    % ------------ 3) number significant (left);
    lt_subplot(3,2,4); hold on;
    title('number sig cases (negative)');
    
    lt_plot_histogram(NsignLeft_ShuffMult, xcenters);
    ydat = sum(rhoPallDAT<0.05 & rhoallDAT<0);
    line([ydat ydat], ylim, 'Color', 'r');
    
    p = (1+sum(NsignLeft_ShuffMult>=ydat))/length(NsignLeft_ShuffMult);
    lt_plot_pvalue(p, 'vs shuff', 1);
    
end



%% %%%%%%%%%%%%%%%%%%%%%%% [DOES FF CORR ACCOUNT FOR CONTEXT DIFF?]
%% ============ PLOT DISTRIBUTION OF ALL MEAN PITCH AND FRATE DIFF
% --------------------- ACTUAL DATA
pitch_as_predictor = 1; % if 1, then pitch on x axis...
doshuff = 0;
[AllBranch_IntDiff, AllBranch_SlopeDiff, AllBranch_SlopeOverall, ...
    AllBranch_IntCoeff_MedAbs, AllBranch_SlopeCoeff_MedAbs, AllBranch_SlopeOverallCoeff_MedAbs, ...
    AllBranch_birdnum, AllBranch_branchnum, AllBranch_neurnum] = ...
    lt_neural_v2_CTXT_Acoustic_CorrSub1(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, pitch_as_predictor, ...
    SummaryStruct, doshuff);

%% 
% AllBranch_IntDiff = AllBranch_IntDiff_DAT;
% AllBranch_SlopeDiff = AllBranch_SlopeDiff_DAT;
% AllBranch_SlopeOverall = AllBranch_SlopeOverall_DAT;
% % AllBranch_birdnum = AllBranch_birdnum_DAT;
% % AllBranch_branchnum = AllBranch_branchnum_DAT;
% % AllBranch_neurnum = AllBranch_neurnum_DAT;


%% ================ EXTRACT DECODE STATS
AllBranch_DecodeP = nan(size(AllBranch_birdnum));
AllBranch_Nctxt = nan(size(AllBranch_birdnum));

for i=1:numbirds
    for bb = 1:maxbranch
        for nn = 1:maxneur
            
            inds = AllBranch_birdnum==i & AllBranch_branchnum==bb & AllBranch_neurnum==nn;
            if ~any(inds)
                continue
            else
                assert(sum(inds)==1, 'asfasd')
            end
                
            % ============== extract decode stats
            DAT = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT;
            indneur = DAT.IndsOriginal_neuron==nn;
            assert(sum(indneur)==1,'asdfas');
            
            decode_p = DAT.PREMOTORDECODE_pval(indneur);
            
            AllBranch_DecodeP(find(inds)) = decode_p;

            % ============ how manyc ontext?
            nctxt = length(DAT.frdat(indneur).classnum);
            AllBranch_Nctxt(find(inds)) = nctxt;
        end
    end
end
assert(~any(isnan(AllBranch_DecodeP)));



%% ============ PLOT PROPORTION CASES SIGNIFICANT CLASS/SLOPE DIFF
lt_figure; hold on;
count = 1;
decodelist = [1 2]; % iterate over diff contingencies on decode performance
plotfraclist = [0 1]; % iterate over plotting either fractions or freq.

for plotfrac = plotfraclist
    % plotfrac = 0; % then instead of actual number cases ,...
    
    for onlydecode = decodelist
        % onlydecode = 1;
        % 0: don't care
        % 1: only cases with significant decode
        % 2: only cases with NON significant decode
        lt_subplot(2,2,count); hold on;
        if plotfrac==1
            ylabel('fraction syls');
        else
            ylabel('n syls');
        end
        
        xlabel('int_only - slope_only - both - neither(slope_main) - neither(no slope_main)');
        for i=1:numbirds
            
            if onlydecode==0
                inds = AllBranch_birdnum==i;
                title('all data (dont care about decode)');
            elseif onlydecode==1
                inds =  AllBranch_birdnum==i & AllBranch_DecodeP<0.05;
                title('only if significant decode');
            elseif onlydecode==2
                inds =  AllBranch_birdnum==i & AllBranch_DecodeP>=0.05;
                title('only if NONsig decode');
            end
            line([i+0.5 i+0.5], ylim);
            
            if ~any(inds)
                continue
            end
            
            intdiff = AllBranch_IntDiff(inds);
            slopediff = AllBranch_SlopeDiff(inds);
            slopeoverall = AllBranch_SlopeOverall(inds);
            
            Y = [];
            % --- only int
            y = sum(intdiff==1 & slopediff==0);
            Y = [Y y];
            
            % --- only slope
            y = sum(intdiff==0 & slopediff==1);
            Y = [Y y];
            
            % --- both
            y = sum(intdiff==1 & slopediff==1);
            Y = [Y y];
            
            % --- none (with overall slope effect)
            y = sum(intdiff==0 & slopediff==0 & slopeoverall==1);
            Y = [Y y];
            
            % --- none (WITHOUT overall slope effect)
            y = sum(intdiff==0 & slopediff==0 & slopeoverall==0);
            Y = [Y y];
            
            % --- nromalize to get fraction
            if plotfrac==1
                Y = Y./length(intdiff);
            end
            
            % --- get x locations
            X = i+[-0.25 -0.1 0.05 0.2 0.35];
            
            % ---------------- SIG EFFECTS
            lt_plot_bar(X(1:3), Y(1:3), {'BarWidth', 0.6, 'Color', 'r'});
            
            % ---------------- NOT SIG (with slope efefct
            lt_plot_bar(X, [nan nan nan Y(4) nan], {'BarWidth', 0.6});
            
            % --------------- NOT SIG (no slope effect)
            lt_plot_bar(X, [nan nan nan nan Y(5)], {'BarWidth', 0.6, 'Color', 'none'});
            
            
        end
        
        % ###################### across entire dataset
            if onlydecode==0
                inds = ones(size(AllBranch_birdnum));
            elseif onlydecode==1
                inds =  AllBranch_DecodeP<0.05;
            elseif onlydecode==2
                inds =  AllBranch_DecodeP>=0.05;
            end
            line([i+0.5 i+0.5], ylim);
            
            intdiff = AllBranch_IntDiff(inds);
            slopediff = AllBranch_SlopeDiff(inds);
            slopeoverall = AllBranch_SlopeOverall(inds);
            
            Y = [];
            % --- only int
            y = sum(intdiff==1 & slopediff==0);
            Y = [Y y];
            
            % --- only slope
            y = sum(intdiff==0 & slopediff==1);
            Y = [Y y];
            
            % --- both
            y = sum(intdiff==1 & slopediff==1);
            Y = [Y y];
            
            % --- none (with overall slope effect)
            y = sum(intdiff==0 & slopediff==0 & slopeoverall==1);
            Y = [Y y];
            
            % --- none (WITHOUT overall slope effect)
            y = sum(intdiff==0 & slopediff==0 & slopeoverall==0);
            Y = [Y y];
            
            % --- nromalize to get fraction
            if plotfrac==1
                Y = Y./length(intdiff);
            end
            
            % --- get x locations
            X = numbirds+1+[-0.25 -0.1 0.05 0.2 0.35];
            
            % ---------------- SIG EFFECTS
            lt_plot_bar(X(1:3), Y(1:3), {'BarWidth', 0.6, 'Color', 'r'});
            
            % ---------------- NOT SIG (with slope efefct
            lt_plot_bar(X, [nan nan nan Y(4) nan], {'BarWidth', 0.6});
            
            % --------------- NOT SIG (no slope effect)
            lt_plot_bar(X, [nan nan nan nan Y(5)], {'BarWidth', 0.6, 'Color', 'none'});
        
        
        % ===================== formatting
        set(gca, 'XTick', 1:numbirds+1, 'XTickLabel', [{SummaryStruct.birds.birdname} 'ALL']);
        rotateXLabels(gca, 90);
        count = count+1;
    end
end

%% =============== [DAT VS. SHUFFLE] EXTRACT;
close all;
% --------------------- PERFORM AND COLLECT MULTIPLE SHUFFLES
pitch_as_predictor = 1; % if 1, then pitch on x axis...
doshuff = 1;
Nshuff = 1000;

AllBranch_SHUFFSTRUCT = struct;
for n=1:Nshuff
    disp(' ------------------------------- ');
    disp(['SHUFFNUM: ' num2str(n)]);
    
    [intdiff_shuff, slopediff_shuff, slopeoverall_shuff, ...
        IntCoeff_MedAbs, SlopeCoeff_MedAbs, SlopeOverallCoeff_MedAbs, ...
        AllBranch_birdnum, AllBranch_branchnum, AllBranch_neurnum] = ...
        lt_neural_v2_CTXT_Acoustic_CorrSub1(AllBirdnum, AllBranchnum, ...
        AllNeurnum, AllClassnum, AllFRmeans, AllPitch, pitch_as_predictor, ...
        SummaryStruct, doshuff);
    
    
    AllBranch_SHUFFSTRUCT.cycle(n).IntDiff = intdiff_shuff;
    AllBranch_SHUFFSTRUCT.cycle(n).SlopeDiff = slopediff_shuff;
    AllBranch_SHUFFSTRUCT.cycle(n).SlopeOverall = slopeoverall_shuff;
    
    AllBranch_SHUFFSTRUCT.cycle(n).IntCoeff_MedAbs = IntCoeff_MedAbs;
    AllBranch_SHUFFSTRUCT.cycle(n).SlopeCoeff_MedAbs = SlopeCoeff_MedAbs;
    AllBranch_SHUFFSTRUCT.cycle(n).SlopeOverallCoeff_MedAbs = SlopeOverallCoeff_MedAbs;
    
end


% ============== SAVE (since takes a long time to do shuffles)
savefname = [savedir '/' analyfname '/Acoustic_Corr/AllBranch_SHUFFSTRUCT'];
if exist([savedir '/' analyfname '/Acoustic_Corr'], 'dir') == 0
    mkdir([savedir '/' analyfname '/Acoustic_Corr']);
end
save(savefname, 'AllBranch_SHUFFSTRUCT');

% ============== SAVE ACTUAL DAT
tmpnames = whos('AllBranch*');
for j=1:length(tmpnames)
savefname = [savedir '/' analyfname '/Acoustic_Corr/' tmpnames(j).name];
save(savefname, tmpnames(j).name);
end

%% ============= [DAT VS> SHUFFLE] SIGNIFICANCE
% FOR EACH CASE, ask how unlikely its effect size is given its own shuffle
% distribution
% ---- which coefficient to care about?
ydat = AllBranch_SlopeCoeff_MedAbs;
yshuff = [AllBranch_SHUFFSTRUCT.cycle.SlopeCoeff_MedAbs];

% ===========================================================
shuffmat = yshuff;
datmat = repmat(ydat, 1, size(shuffmat,2));

% ------ get vals
pvalmat = (sum(shuffmat>=datmat, 2))./size(shuffmat,2);

ind_sigcases = pvalmat<0.01;
disp([num2str(sum(ind_sigcases)) ' sig cases out of ' num2str(length(ind_sigcases))]);

% ------ zscore all data relative to shuffle distribution
lt_figure; hold on;
lt_subplot(3,2,1); hold on;
title('singal case shuff distr (rand choice)');
casethis = randi(size(yshuff,1),1);
lt_plot_histogram(yshuff(casethis,:));

lt_subplot(3,2,2); hold on;
title('all case zscore rel own shuffles');
y_zscore = (ydat - mean(yshuff,2))./std(yshuff,0,2);

lt_plot_histogram(y_zscore);
lt_plot_zeroline_vert;

lt_subplot(3,2,3); hold on;
title('all cases, percentiles rel own shuffle');
xcenters = 0:0.05:1;
[Ybinned] = lt_plot_histogram(pvalmat, xcenters);

lt_subplot(3,2,4); hold on;
title('all cases, percentiles rel own shuffle');
xlabel('prct');
ylabel('frac cases equal to or below this percent')

x = xcenters;
binsize = xcenters(2)-xcenters(1);
x = x+binsize/2;
x(end) = 1;
y = cumsum(Ybinned)./sum(Ybinned);
plot(x,y, '-k');
line([0 1], [0 1], 'Color', 'b', 'LineStyle', '--');

%% ============= [DAT VS. SHUFFLE] COMPARE DISTRIBUTIONS
lt_figure; hold on;
shuffcycle = 1; % which one ot plot;
onlysigdecode = 0; % then filters based on decode [0 1 or 2]
NctxtToPlot = [1:10]; % can't be empty 

if onlysigdecode==1
    inds = AllBranch_DecodeP<0.05 & ismember(AllBranch_Nctxt, NctxtToPlot);
elseif onlysigdecode==0
    inds = logical(ones(size(AllBranch_DecodeP))) & ismember(AllBranch_Nctxt, NctxtToPlot);
elseif onlysigdecode==2
    inds = AllBranch_DecodeP>=0.05 & ismember(AllBranch_Nctxt, NctxtToPlot);
end

% ---- NCASES, slope interaction
lt_subplot(4,2,1); hold on;
title('ncases slope interaction');
y = sum(AllBranch_SlopeDiff(inds));
tmp = [AllBranch_SHUFFSTRUCT.cycle.SlopeDiff];
yshuff = sum(tmp(inds,:),1);

lt_plot_histogram(yshuff);
line([y, y], ylim, 'Color', 'r');
lt_plot_annotation(1, ['N=' num2str(sum(inds))], 'b');
% pval
p = (1+sum(yshuff>=y))./length(yshuff);
lt_plot_pvalue(p, 'vs. shuff', 1);


% ---- [SINGLE SHUFF] median absolute(slope-interaction effect)
lt_subplot(4,2,2); hold on;
title('[SINGLE SHUFF]median(across classes) abs(slope-context interact)');
xlabel('bk=shuff');
y = AllBranch_SlopeCoeff_MedAbs(inds);
yshuff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeCoeff_MedAbs(inds);
[~, xcenters] = lt_plot_histogram(y, '', 1, 1, '', 1, 'r');
lt_plot_histogram(yshuff, xcenters, 1, 1, '', 1, 'k');

p = signrank(y, yshuff);
lt_plot_pvalue(p, 'srank', 1);

% ---- [MULT SHUFF, each with median] median absolute(slope-interaction effect)
lt_subplot(4,2,3); hold on;
title('[MULT SHUFF]median(across classes) abs(slope-context interact)');

y = median(AllBranch_SlopeCoeff_MedAbs(inds));
tmp = [AllBranch_SHUFFSTRUCT.cycle.SlopeCoeff_MedAbs];
yshuff = median(tmp(inds,:), 1);
lt_plot_histogram(yshuff, '', 1, 1, '', 1, 'k');
line([y y], ylim, 'Color', 'r');
% -pvalue
p = (1+sum(yshuff>=y))./length(yshuff);
lt_plot_pvalue(p, 'vs. shuff', 1);
lt_plot_zeroline_vert;


% ---- [SINGLE SHUFF] overall slopes
lt_subplot(4,2,4); hold on;
title('[SINGLE SHUFF] median(across classes) abs(overall slope)');
xlabel('bk=shuff');
y = AllBranch_SlopeOverallCoeff_MedAbs(inds);
yshuff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeOverallCoeff_MedAbs(inds);
[~, xcenters] = lt_plot_histogram(y, '', 1, 1, '', 1, 'r');
lt_plot_histogram(yshuff, xcenters, 1, 1, '', 1, 'k');

p = signrank(y, yshuff);
lt_plot_pvalue(p, 'srank', 1);


% ---- [SINGLE SHUFF] intercept x class interaction
lt_subplot(4,2,5); hold on;
title('[SINGLE SHUFF] median(across classes) abs(intercept-context interaction)');
xlabel('bk=shuff');
y = AllBranch_IntCoeff_MedAbs(inds);
yshuff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).IntCoeff_MedAbs(inds);
[~, xcenters] = lt_plot_histogram(y, '', 1, 1, '', 1, 'r');
lt_plot_histogram(yshuff, xcenters, 1, 1, '', 1, 'k');

p = signrank(y, yshuff);
lt_plot_pvalue(p, 'srank', 1);



%% ============= [DAT VS. SHUFFLE] PLOT
%% [PLOT] OVERLAY DATA WITH ONE INSTANCE OF SHUFFLE
lt_figure; hold on;
count = 1;
decodelist = [1 2]; % iterate over diff contingencies on decode performance
plotfrac = 1;
shuffcycle = 1;
for useshuff = [0 1]
    for onlydecode = decodelist
        % onlydecode = 1;
        % 0: don't care
        % 1: only cases with significant decode
        % 2: only cases with NON significant decode
        lt_subplot(2,2,count); hold on;
        if plotfrac==1
            ylabel('fraction syls');
        else
            ylabel('n syls');
        end
        
        xlabel('int_only - slope_only - both - neither(slope_main) - neither(no slope_main)');
        for i=1:numbirds
            
            if onlydecode==0
                inds = AllBranch_birdnum==i;
                title('all data (dont care about decode)');
            elseif onlydecode==1
                inds =  AllBranch_birdnum==i & AllBranch_DecodeP<0.05;
                title('only if significant decode');
            elseif onlydecode==2
                inds =  AllBranch_birdnum==i & AllBranch_DecodeP>=0.05;
                title('only if NONsig decode');
            end
            line([i+0.5 i+0.5], ylim);
            
            if ~any(inds)
                continue
            end
            
            if useshuff==0
                % then use dat
                intdiff = AllBranch_IntDiff(inds);
                slopediff = AllBranch_SlopeDiff(inds);
                slopeoverall = AllBranch_SlopeOverall(inds);
            else
                intdiff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).IntDiff(inds);
                slopediff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeDiff(inds);
                slopeoverall = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeOverall(inds);
            end
            
            Y = [];
            % --- only int
            y = sum(intdiff==1 & slopediff==0);
            Y = [Y y];
            
            % --- only slope
            y = sum(intdiff==0 & slopediff==1);
            Y = [Y y];
            
            % --- both
            y = sum(intdiff==1 & slopediff==1);
            Y = [Y y];
            
            % --- none (with overall slope effect)
            y = sum(intdiff==0 & slopediff==0 & slopeoverall==1);
            Y = [Y y];
            
            % --- none (WITHOUT overall slope effect)
            y = sum(intdiff==0 & slopediff==0 & slopeoverall==0);
            Y = [Y y];
            
            % --- nromalize to get fraction
            if plotfrac==1
                Y = Y./length(intdiff);
            end
            
            % --- get x locations
            X = i+[-0.25 -0.1 0.05 0.2 0.35];
            
            % ---------------- SIG EFFECTS
            lt_plot_bar(X(1:3), Y(1:3), {'BarWidth', 0.6, 'Color', 'r'});
            
            % ---------------- NOT SIG (with slope efefct
            lt_plot_bar(X, [nan nan nan Y(4) nan], {'BarWidth', 0.6});
            
            % --------------- NOT SIG (no slope effect)
            lt_plot_bar(X, [nan nan nan nan Y(5)], {'BarWidth', 0.6, 'Color', 'none'});
            
            
        end
        
        % ###################### across entire dataset
        if onlydecode==0
            inds = ones(size(AllBranch_birdnum));
        elseif onlydecode==1
            inds =  AllBranch_DecodeP<0.05;
        elseif onlydecode==2
            inds =  AllBranch_DecodeP>=0.05;
        end
        line([i+0.5 i+0.5], ylim);
        
            if useshuff==0
                % then use dat
                intdiff = AllBranch_IntDiff(inds);
                slopediff = AllBranch_SlopeDiff(inds);
                slopeoverall = AllBranch_SlopeOverall(inds);
            else
                intdiff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).IntDiff(inds);
                slopediff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeDiff(inds);
                slopeoverall = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeOverall(inds);
            end
        
        Y = [];
        % --- only int
        y = sum(intdiff==1 & slopediff==0);
        Y = [Y y];
        
        % --- only slope
        y = sum(intdiff==0 & slopediff==1);
        Y = [Y y];
        
        % --- both
        y = sum(intdiff==1 & slopediff==1);
        Y = [Y y];
        
        % --- none (with overall slope effect)
        y = sum(intdiff==0 & slopediff==0 & slopeoverall==1);
        Y = [Y y];
        
        % --- none (WITHOUT overall slope effect)
        y = sum(intdiff==0 & slopediff==0 & slopeoverall==0);
        Y = [Y y];
        
        % --- nromalize to get fraction
        if plotfrac==1
            Y = Y./length(intdiff);
        end
        
        % --- get x locations
        X = numbirds+1+[-0.25 -0.1 0.05 0.2 0.35];
        
        % ---------------- SIG EFFECTS
        lt_plot_bar(X(1:3), Y(1:3), {'BarWidth', 0.6, 'Color', 'r'});
        
        % ---------------- NOT SIG (with slope efefct
        lt_plot_bar(X, [nan nan nan Y(4) nan], {'BarWidth', 0.6});
        
        % --------------- NOT SIG (no slope effect)
        lt_plot_bar(X, [nan nan nan nan Y(5)], {'BarWidth', 0.6, 'Color', 'none'});
        
        
        % ===================== formatting
        set(gca, 'XTick', 1:numbirds+1, 'XTickLabel', [{SummaryStruct.birds.birdname} 'ALL']);
        rotateXLabels(gca, 90);
        count = count+1;
        if useshuff==1
            lt_plot_annotation(1, 'SHUFF', 'b');
        end
    end
end

%% = PLOT "RASTER" OF ALL CASES, INDICATING MAGNITUDE AND SIGNIFICANCE
effecttosortby = 2;

plotshuffle = 1; % have to have done shuffle already
sc = 1; % choose which one.

% ================
lt_figure; hold on;
ylabel('case #');
xlabel('Int*ctxt -- Slope*ctxt -- SlopeOverall');

% --- ciollect data
if plotshuffle==0
Y = [AllBranch_IntCoeff_MedAbs AllBranch_SlopeCoeff_MedAbs AllBranch_SlopeOverallCoeff_MedAbs];
Ysig = [AllBranch_IntDiff AllBranch_SlopeDiff AllBranch_SlopeOverall];
Ydecode = AllBranch_DecodeP<0.05;
title('[DAT]median abs(cofficeint)');
elseif plotshuffle==1
Y = [AllBranch_SHUFFSTRUCT.cycle(sc).IntCoeff_MedAbs AllBranch_SHUFFSTRUCT.cycle(sc).SlopeCoeff_MedAbs ...
    AllBranch_SHUFFSTRUCT.cycle(sc).SlopeOverallCoeff_MedAbs];
Ysig = [AllBranch_SHUFFSTRUCT.cycle(sc).IntDiff AllBranch_SHUFFSTRUCT.cycle(sc).SlopeDiff ...
    AllBranch_SHUFFSTRUCT.cycle(sc).SlopeOverall];
Ydecode = AllBranch_DecodeP<0.05;
title('[SHUFF]median abs(cofficeint)');
end
% ---- sort data, increasing effect *(choose effect)
[~, inds] = sort(Y(:,effecttosortby));
Y = Y(inds,:);
Ysig = Ysig(inds, :);
Ydecode = Ydecode(inds,:);


% ---- plot
x = 1:size(Y,2);
imagesc(Y, [0 2]);
colorbar;
colormap('gray');

% ---- mark those that are significnat
for j=1:size(Ysig,2)
   plot(j-0.3, find(Ysig(:,j)), 'or'); 
end
    
% ---- mark those with significnat decode (regression)
plot(0.3, find(Ydecode), 'ob');    

% lt_figure; hold on;
% for j=1:size(Y,2)
%     lt_plot_stem3(j*ones(size(Y,1),1), 1:size(Y,1)', Y(:, j), 'k', 1);
% end


%% ========= [PLOT] plot ff vs. fr (overlay in diff color for all classes)
% OVERLAY OUTCOME OF CONJUNCTIVE CODING

for i=1:numbirds
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    birdname = SummaryStruct.birds(i).birdname;
    
    for bb = 1:maxbranch
        
        for nn = 1:maxneur
            
            inds = AllBirdnum==i & AllBranchnum==bb & AllNeurnum==nn;
            if ~any(inds)
                continue
            end
            
            branchname = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            % ----------------
            classthis = AllClassnum{inds};
            fratemeanthis = sqrt(AllFRmeans{inds});
            pitchthis = AllPitch{inds};
            
            if all(isnan(pitchthis))
                continue
            end
            
            % ========= plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-' branchname{1} num2str(bb) '-n' num2str(nn)]);
            if pitch_as_predictor==1
                xlabel('pitch');
                ylabel('mean frate (sqrt)');
            else
                ylabel('pitch');
                xlabel('mean frate (sqrt)');
            end
            assert(sum(inds)==1,'asdfsd');
            
            % ---------------------------
            pcols = lt_make_plot_colors(max(classthis), 0,0);
            for cc=1:max(classthis)
                indtmp = classthis==cc;
                
                if pitch_as_predictor==1
                    x = pitchthis(indtmp);
                    y = fratemeanthis(indtmp);
                else
                    % then frate on x axis
                    x = fratemeanthis(indtmp);
                    y = pitchthis(indtmp);
                end
                
                plot(x,y, 'x', 'Color', pcols{cc});
                
                % -------- plot mean pitch and frate for this class
                xmean = mean(x);
                xsem = lt_sem(x);
                ymean = mean(y);
                ysem = lt_sem(y);
                
                %                lt_plot(xmean, ymean, {'Errors', ysem, 'Color', pcols{cc}});
                lt_plot(xmean, ymean, {'Color', pcols{cc}});
                line(xmean+[-xsem xsem], [ymean ymean], 'Color', pcols{cc});
                line([xmean xmean], ymean+[-ysem ysem], 'Color', pcols{cc});
            end
            
            % ============ analysis outocme
            indbranch = AllBranch_birdnum==i & AllBranch_branchnum==bb & ...
                AllBranch_neurnum==nn;
            assert(sum(indbranch)==1, 'asfdas');
            intdiff = AllBranch_IntDiff(indbranch);
            slopediff = AllBranch_SlopeDiff(indbranch);
            
            lt_plot_annotation(1, ['intdiff' num2str(intdiff) '; slope' num2str(slopediff)], 'r');
            
        end
    end
    
    pause;
    close all;
end

end


%% FUNCTIONS

function [rhoall, rhoCIall, rhoPall, Nsampall] = fn_getallRho(AllBirdnum, AllBranchnum, ...
    AllNeurnum, AllClassnum, AllFRmeans, AllPitch, doshuff, poolOverCtxt)

if ~exist('poolOverCtxt', 'var')
    % if 1, then ignore context (e.g. Sober analysis)
    poolOverCtxt=0;
end
%%
numbirds = max(AllBirdnum);
maxbranch = max(AllBranchnum);
maxneur = max(AllNeurnum);

% =====================
rhoall =[];
rhoCIall = [];
Nsampall =[];
rhoPall = [];

for i=1:numbirds
    %     birdname = SummaryStruct.birds(i).birdname;
    for bb = 1:maxbranch
        for nn = 1:maxneur
            
            inds = AllBirdnum==i & AllBranchnum==bb & AllNeurnum==nn;
            if ~any(inds)
                continue
            end
            assert(sum(inds)==1,'sdfasdf');
            
            % ----------------
            classthis = AllClassnum{inds};
            fratemeanthis = sqrt(AllFRmeans{inds});
            pitchthis = AllPitch{inds};
            
            if all(isnan(pitchthis))
                continue
            end
            
            
            % ===================== SEPARATING BY CONTEXT
            if poolOverCtxt==0
                % ---------------------------
                %             pcols = lt_make_plot_colors(max(classthis), 0,0);
                for cc=1:max(classthis)
                    indtmp = classthis==cc;
                    
                    x = fratemeanthis(indtmp);
                    y = pitchthis(indtmp);
                    
                    % ================== shuffle trials (within individual
                    % datasets)?
                    if doshuff==1
                        % break link between neural and beahvior
                        y = y(randperm(length(y)));
                    end
                    
                    % ================ get corr with CI
                    %                 rho = corr(x',y');
                    [rho, p, rhoL, rhoU] = corrcoef(x', y');
                    rho = rho(1,2);
                    rhoCI = [rhoL(1,2) rhoU(1,2)];
                    p = p(1,2);
                    
                    % ======= save
                    rhoall =[rhoall; rho];
                    rhoCIall = [rhoCIall; rhoCI];
                    Nsampall = [Nsampall; length(x)];
                    rhoPall = [rhoPall; p];
                end
            else
                % ================== POOL OVER CONTEXT
                    x = fratemeanthis;
                    y = pitchthis;
                    
                    % ================== shuffle trials (within individual
                    % datasets)?
                    if doshuff==1
                        % break link between neural and beahvior
                        y = y(randperm(length(y)));
                    end
                    
                    % ================ get corr with CI
                    %                 rho = corr(x',y');
                    [rho, p, rhoL, rhoU] = corrcoef(x', y');
                    rho = rho(1,2);
                    rhoCI = [rhoL(1,2) rhoU(1,2)];
                    p = p(1,2);
                    
                    % ======= save
                    rhoall =[rhoall; rho];
                    rhoCIall = [rhoCIall; rhoCI];
                    Nsampall = [Nsampall; length(x)];
                    rhoPall = [rhoPall; p];

            end
        end
    end
end

end





















