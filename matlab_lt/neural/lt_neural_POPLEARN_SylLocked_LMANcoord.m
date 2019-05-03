function DATSTRUCT_POP = lt_neural_POPLEARN_SylLocked_LMANcoord(DATSTRUCT_POP, DATSTRUCT_SPK, ...
    DATSTRUCT_LFP, PARAMS, OUTSTRUCT_XCOV, SwitchStruct, dirstruct, OUTSTRUCT, ...
    bregionthis, plotbase)
%% ======== LT 4/5/19 - Hamish-style LMAN coordination, how relate to things:
% 1) from trial to trial relate to LMAN-RA xcov?
% 2) correlate with pitch?
% 3) pitch correlation relates to directed song bias?

%% ======== FOR EACH CASE, GET CORRELATION BETWEEN lman ENSEMBLE AND PITCH

SlopeCI_all = cell(size(DATSTRUCT_POP.bnum,1),1);

for i=1:length(DATSTRUCT_POP.bnum)
    
    indsthis = i;
    
    bnum = DATSTRUCT_POP.bnum(indsthis);
    enum = DATSTRUCT_POP.enum(indsthis);
    sw = DATSTRUCT_POP.switch(indsthis);
    mm = DATSTRUCT_POP.motifnum (indsthis);
    
    % ==== correaltion with pitch
    rhoEnsemble = DATSTRUCT_POP.fr_RhoPairwise{indsthis};
    indsbase = DATSTRUCT_POP.indsbase_epoch{indsthis};
    indsWN = DATSTRUCT_POP.indsWN_epoch{indsthis};
    if plotbase==1
        indstouse = indsbase;
    else
        indstouse = indsWN;
    end
        
    % -- get FF
    indsLFP = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw ...
        & OUTSTRUCT.motifnum==mm);
    indsLFP = indsLFP(1);
    
%     tvals = OUTSTRUCT.tvals{indsLFP};
    ffvals = OUTSTRUCT.ffvals{indsLFP};
%     motifname = OUTSTRUCT.motifname{indsLFP};
    
    assert(length(rhoEnsemble)==length(ffvals), 'why no match? not sure');
    
    
    % ======== SKIP IF NO DEFINED PITCH
    if all(isnan(ffvals))
        continue
    end
    
    
    f = ffvals(indstouse);
    rho = rhoEnsemble(indstouse);
    [b,bint,r,rint,stats,SummaryStats] = ...
        lt_regress(f, rho, 0, 0, 1, 1, 'k');
   
    slope_CI = [SummaryStats.slopeCI(1) SummaryStats.slope SummaryStats.slopeCI(2)];    
    SlopeCI_all{i} = slope_CI;
    
end

DATSTRUCT_POP.Slope_Pitch_vs_LMANcorr_CI = SlopeCI_all;


%% ======== correlate with pitch?

[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_POP.bnum, DATSTRUCT_POP.enum, ...
    DATSTRUCT_POP.switch, DATSTRUCT_POP.motifnum});

figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

if strcmp(bregionthis, 'LMAN')
    pcol = [0.3 0.7 0.3];
elseif strcmp(bregionthis, 'RA')
    pcol = [0.8 0.2 0.2];
end
    
for i=1:length(indsgrpU)
   indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, bregionthis));
   
   if isempty(indsthis)
       continue
   end
   
   assert(length(indsthis)==1);
   
   bnum = DATSTRUCT_POP.bnum(indsthis);
   enum = DATSTRUCT_POP.enum(indsthis);
   sw = DATSTRUCT_POP.switch(indsthis);
   mm = DATSTRUCT_POP.motifnum(indsthis);
   slopethis = DATSTRUCT_POP.Slope_Pitch_vs_LMANcorr_CI{indsthis};
   
   % ==== correaltion with pitch
    rhoEnsemble = DATSTRUCT_POP.fr_RhoPairwise{indsthis};
    indsbase = DATSTRUCT_POP.indsbase_epoch{indsthis};
    indsWN = DATSTRUCT_POP.indsWN_epoch{indsthis};
    indstouse = indsbase;
    
    % -- get FF
    indsLFP = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw ...
        & OUTSTRUCT.motifnum==mm);
    indsLFP = indsLFP(1);
    
    tvals = OUTSTRUCT.tvals{indsLFP};
    ffvals = OUTSTRUCT.ffvals{indsLFP};
    motifname = OUTSTRUCT.motifname{indsLFP};
    
    assert(length(rhoEnsemble)==length(ffvals), 'why no match? not sure');
    
    
    % ======== SKIP IF NO DEFINED PITCH
    if all(isnan(ffvals))
        continue
    end

    
    % ================ PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([num2str(bnum) '-' num2str(enum) '-' num2str(sw) '-' motifname], 'Color', pcol);
    xlabel([bregionthis ' ensemble corr']);
    ylabel('ff');
    
    f = ffvals(indstouse);
    rho = rhoEnsemble(indstouse);
    lt_regress(f, rho, 1, 0, 1, 1, 'k');
    
end


%% ================== SUMMARY PLOT OF SLOPES AND CI


%% ================== SUMMARY: PITCH-LMAN corr vs. AFP bias
% ======= 1) COLLECT all AFP biases
AfpBiasAll = nan(size(DATSTRUCT_POP.bnum,1),1);

for i=1:length(DATSTRUCT_POP.bnum)
   
    % -------- find metadat
   bnum = DATSTRUCT_POP.bnum(i);
   bname = SwitchStruct.bird(bnum).birdname;
   enum = DATSTRUCT_POP.enum(i);
   sw = DATSTRUCT_POP.switch(i);
   mm = DATSTRUCT_POP.motifnum(i);
   
    indsLFP = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==sw ...
        & OUTSTRUCT.motifnum==mm);
    indsLFP = indsLFP(1);
    
    mID = OUTSTRUCT.motifID_unique(indsLFP);
    
    % ---------
    inddir = find(strcmp({dirstruct.bird.birdname}, bname));
    afpbias = dirstruct.bird(inddir).DAT.motifID(mID).means.afpbias;
    
    
    % =============== PLOT
    AfpBiasAll(i) = afpbias;
    
end

DATSTRUCT_POP.AFPbias = AfpBiasAll;
    
%% ============= RELATIONSHIP BETWEEN AFP BIAS AND NEURAL?

[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_POP.bnum, DATSTRUCT_POP.enum, ...
    DATSTRUCT_POP.switch});


figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

if strcmp(bregionthis, 'LMAN')
    pcol = [0.3 0.7 0.3];
elseif strcmp(bregionthis, 'RA')
    pcol = [0.8 0.2 0.2];
end
    
for i=1:length(indsgrpU)
   indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, bregionthis));
   
   if isempty(indsthis)
       continue
   end
   
   
   bnum = unique(DATSTRUCT_POP.bnum(indsthis));
   enum = unique(DATSTRUCT_POP.enum(indsthis));
   sw = unique(DATSTRUCT_POP.switch(indsthis));
   
    % =========== stats
    slope = DATSTRUCT_POP.Slope_Pitch_vs_LMANcorr_CI(indsthis);
    afpbias = DATSTRUCT_POP.AFPbias(indsthis);
    
    
    % ================ PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([num2str(bnum) '-' num2str(enum) '-' num2str(sw)], 'Color', pcol);
    xlabel('slope (PitchVsLMANcorr, hz/rho)');
    ylabel('afpbias (hz)');
    hsplots = [hsplots hsplot];
    
    % ------------ 1) plot with error bars
    for j=1:length(slope)
        if isempty(slope{j}) 
            continue
        end
        plot(slope{j}(2), afpbias(j), 'ok');
        line([slope{j}(1) slope{j}(3)], [afpbias(j) afpbias(j)], 'Color', 'k');
    end
    
    
%     lt_regress(f, rho, 1, 0, 1, 1, 'k');
    lt_plot_zeroline;
    lt_plot_zeroline_vert
end

linkaxes(hsplots, 'xy');


% ############################### PLOT ALL ON SAME PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('all same plot', 'Color', pcol);
    xlabel('slope (PitchVsLMANcorr, hz/rho)');
    ylabel('afpbias (hz)');
    hsplots = [hsplots hsplot];
    pcols = lt_make_plot_colors(length(indsgrpU), 0,0);
    
    for i=1:length(indsgrpU)
   indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, bregionthis));
   
   if isempty(indsthis)
       continue
   end
   
   
   bnum = unique(DATSTRUCT_POP.bnum(indsthis));
   enum = unique(DATSTRUCT_POP.enum(indsthis));
   sw = unique(DATSTRUCT_POP.switch(indsthis));
   
    % =========== stats
    slope = DATSTRUCT_POP.Slope_Pitch_vs_LMANcorr_CI(indsthis);
    afpbias = DATSTRUCT_POP.AFPbias(indsthis);
    
    
    % ================ PLOT
    % ------------ 1) plot with error bars
    for j=1:length(slope)
        if isempty(slope{j}) 
            continue
        end
%         lt_plot(slope{j}(2), afpbias(j), {'Color', pcols{i}});
        plot(slope{j}(2), afpbias(j), 'o', 'Color', pcols{i});
%         line([slope{j}(1) slope{j}(3)], [afpbias(j) afpbias(j)], 'Color', 'k');
        line([slope{j}(1) slope{j}(3)], [afpbias(j) afpbias(j)], 'Color', pcols{i});
    end
    
    
%     lt_regress(f, rho, 1, 0, 1, 1, 'k');
    lt_plot_zeroline;
    lt_plot_zeroline_vert
end

linkaxes(hsplots, 'xy');


%% ========================== 


