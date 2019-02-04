%% ===================== [PLOT] EXTRACT MOTIFS (FROM LABEL);
% 
% PARAMS.MotifsToGet = {...
%     'pu69wh78', {'(a)abhh', '(j)jbhh'}, ...
%     'wh72pk12', {}, ...
%     };
PARAMS.MotifsToGet = {...
    'pu69wh78', {'aa(b)hh', 'jj(b)hh'}, ...
    'wh72pk12', {}, ...
    };

PARAMS.motif_predur = 0.3;
PARAMS.motif_postdur = 0.2; %
PARAMS.usemotifoffset = 1; % if 0, then from onset of token.

DATSTRUCT = struct;
DATSTRUCT.Spktimes = {};
DATSTRUCT.Onstimes= {};
DATSTRUCT.Offtimes= {};
DATSTRUCT.unitnum = [];
DATSTRUCT.exptnum = [];
DATSTRUCT.songrendnum = [];
DATSTRUCT.motifnum = [];
DATSTRUCT.motifname = {};
DATSTRUCT.motifregexp = {};
DATSTRUCT.birdname = {};
DATSTRUCT.bostype = [];
DATSTRUCT.minmax_time = [];
DATSTRUCT.bostype_str = {};
DATSTRUCT.brainregion = {};


% ============= for each unit, extract all cases of each motif
for i=1:length(SummaryBOS.expt)
    birdname = SummaryBOS.expt(i).birdname;
    motifstoget = PARAMS.MotifsToGet{find(strcmp(PARAMS.MotifsToGet, birdname))+1};
    BOSlabels = PARAMS.BOSlabels{strcmp(PARAMS.BOSbirdname, birdname)};
    BOSnames = PARAMS.BOSnames{strcmp(PARAMS.BOSbirdname, birdname)};
    
    % ====== go thru all song renditions
    nrends = length(SummaryBOS.expt(i).DAT_bysongrend.BOStype);
    
    for nn=1:nrends
        spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset(nn,:);
        bostype = SummaryBOS.expt(i).DAT_bysongrend.BOStype(nn);
        bostype_str = BOSnames{bostype};
        labels = BOSlabels{bostype};
        ons = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets{nn};
        offs = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets{nn};
        
        % ============ go thru all motifs
        for mm=1:length(motifstoget)
            motifthis = motifstoget{mm};
            motiflength = length(motifthis)-2; assert(length(strfind(motifthis, '('))==1, 'need to use parantheses on tokens');
            [startinds, tokeninds, endinds, matchlabs] = lt_neural_QUICK_regexp(labels, motifthis);
            
            % ======= extract data [iterate over all matched labels]
            for jj=1:length(tokeninds)
                
                tt = tokeninds(jj);
                
                datonset = ons(tt)-PARAMS.motif_predur;
                if PARAMS.usemotifoffset ==1
                    datoffset = offs(endinds(jj))+PARAMS.motif_postdur;
                else
                    datoffset = ons(tt)+PARAMS.motif_postdur;
                end
                
                % ==== collect onset and offset times (token onset =0);
                onsthis = ons(startinds(jj):endinds(jj)) - ons(tt);
                offthis = offs(startinds(jj):endinds(jj)) - ons(tt);
                
                % ==== collect spike times
                % ---- once for each chan
                nchans = length(spks);
                for cc=1:nchans
                    spkthis = spks{cc};
                    spkthis = spkthis(spkthis>datonset & spkthis<datoffset);
                    % subtract onset of first syl;
                    spkthis = spkthis - ons(tt);
                    
                    % ============= COLLECT DATA
%                     if isempty(spkthis)
                        DATSTRUCT.Spktimes = [DATSTRUCT.Spktimes; {spkthis}];
%                     else
%                         DATSTRUCT.Spktimes = [DATSTRUCT.Spktimes; spkthis];
%                     end
                    DATSTRUCT.Onstimes= [DATSTRUCT.Onstimes; onsthis];
                    DATSTRUCT.Offtimes= [DATSTRUCT.Offtimes; offthis];
                    DATSTRUCT.unitnum = [DATSTRUCT.unitnum; cc];
                    DATSTRUCT.exptnum = [DATSTRUCT.exptnum; i];
                    DATSTRUCT.songrendnum = [DATSTRUCT.songrendnum; nn];
                    DATSTRUCT.motifnum = [DATSTRUCT.motifnum; mm];
                    DATSTRUCT.motifname = [DATSTRUCT.motifname; matchlabs{jj}];
                    DATSTRUCT.motifregexp = [DATSTRUCT.motifregexp; motifthis];
                    DATSTRUCT.birdname = [DATSTRUCT.birdname; birdname];
                    DATSTRUCT.bostype = [DATSTRUCT.bostype; bostype];
                    DATSTRUCT.bostype_str = [DATSTRUCT.bostype_str; bostype_str];
                    
                    datoffset-datonset-PARAMS.motif_predur;
                    DATSTRUCT.minmax_time = [DATSTRUCT.minmax_time; ...
                        [-PARAMS.motif_predur datoffset-datonset-PARAMS.motif_predur]];
                    
                    DATSTRUCT.brainregion = [DATSTRUCT.brainregion; ...
                        SummaryBOS.expt(i).bregions{cc}];
                end
            end
        end
    end
end

% ========= LIST ALL THINGS EXTRACTED

[inds_out, inds_unique, X_cell] = ...
    lt_tools_grp2idx({DATSTRUCT.birdname DATSTRUCT.exptnum DATSTRUCT.unitnum DATSTRUCT.motifname});
cellfun(@disp, X_cell);


%% ====================== [PLOT]
bird = 'pu69wh78';
motif = 'aabhh';

indsthis = strcmp(DATSTRUCT.birdname, bird) & strcmp(DATSTRUCT.motifname, motif);



%% ========================= [TO DO] LINEAR TIME WARP FUNCTION.


%% ======================== [PLOT] - FOR EACH UNIT, PLOT MEAN MOTIF RESPONSE
close all;

overlayAllMotifs = 1; % one plot per unit, overlaying motifs (from all BOS files)

% --- iterate over experiments and units.
[indsgrp, indsgrp_unique] = lt_tools_grp2idx({DATSTRUCT.exptnum DATSTRUCT.unitnum});

motifslist = unique(DATSTRUCT.motifnum);
bostypelist = unique(DATSTRUCT.bostype);

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



for i=1:length(indsgrp_unique)
    
    pcol = 0.8-0.8*[rand rand rand];
    hsplots = [];
    % =========== go thru all motifs & bostypes
    if overlayAllMotifs==1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         hsplots = [hsplots hsplot];
    end
    for mm = motifslist'
        for bb=1:length(bostypelist)
            bosthis = bostypelist(bb);
            
            indsthis = indsgrp==indsgrp_unique(i) & DATSTRUCT.motifnum==mm & DATSTRUCT.bostype==bosthis;
            if isempty(indsthis)
                continue
            end
            
            pcolmot = 0.8-0.8*[rand rand rand];
            %           bostype = DATSTRUCT.bostype(indsthis);
            rendnum = DATSTRUCT.songrendnum(indsthis);
            spks = DATSTRUCT.Spktimes(indsthis);
            sylonsets = cell2mat(DATSTRUCT.Onstimes(indsthis));
            syloffsets = cell2mat(DATSTRUCT.Offtimes(indsthis));
            
            minmaxtime = DATSTRUCT.minmax_time(indsthis,:);
            
            % ========= currently assuming don't care about type of BOS.
            %             assert(length(unique(bostype))==1);
            
            % === for pplotting.
            birdname = unique(DATSTRUCT.birdname(indsthis));
            exptnum = unique(DATSTRUCT.exptnum(indsthis));
            unitnum = unique(DATSTRUCT.unitnum(indsthis));
            bregion = SummaryBOS.expt(exptnum).bregions{unitnum};
            motif = unique(DATSTRUCT.motifregexp(indsthis));
            bostype = PARAMS.BOSnames{strcmp(PARAMS.BOSbirdname, birdname)}{bosthis};
            
            %% ========== PLOT FOR THIS UNIT, RESPONSES OVER ALL RENDITIONS OF
            % THIS MOTIF
            assert(length(unique(minmaxtime(:,1)))==1, 'assumes all predur from same alignment');
            mintime = minmaxtime(1,1);
            maxtime = min(minmaxtime(:,2));
            
            
            [FRmat, t] = lt_neural_QUICK_Spk2FRmat(spks, mintime, maxtime);
            
            frmean = mean(FRmat,2);
            frsem = lt_sem(FRmat');
            
            
            if overlayAllMotifs==0
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title({[birdname{1} '-enum' num2str(exptnum) '-unit' num2str(unitnum) '[' bregion ']'], ...
                    [motif{1} '-BOS"' bostype '"']});
            else
                lt_plot_text(t(end), frmean(end),  [motif{1} '-BOS"' bostype '"'], pcolmot);
            end

            
            % --- SYLS ONSET OFFSET
            on = median(sylonsets,1);
            off = median(syloffsets,1);
            YLIM = [0 max(frmean)+5];
            
            if overlayAllMotifs==0
                lt_neural_QUICK_PlotSylPatches(on, off, YLIM, 0, pcol);
            else
                lt_neural_QUICK_PlotSylPatches(on, off, YLIM, 1, pcolmot);
            end
            
            
            % --- plot fr mean
             if overlayAllMotifs==0
           shadedErrorBar(t, frmean, frsem, {'Color', pcol},1);
            axis tight;
            lt_plot_zeroline;
             else
           shadedErrorBar(t, frmean, frsem, {'Color', pcolmot},1);
             end
            
        end
    end
    if overlayAllMotifs==1
        title([birdname{1} '-enum' num2str(exptnum) '-unit' num2str(unitnum) '[' bregion ']']);
                axis tight;
            lt_plot_zeroline;
else
    linkaxes(hsplots, 'xy');
    end
end



%% ====================== [ANALYSIS - eXTRACT, EACH UNIT ONE DATAPOINT]
% === COLLECT MEAN RESPONSE TO A GIVEN MOTIF FROM A GIVEN BOS TYPE ACROSS
% ALL UNITS

% --- iterate over experiments and units.
[indsgrp, indsgrp_unique] = lt_tools_grp2idx({DATSTRUCT.exptnum DATSTRUCT.unitnum});

motifslist = unique(DATSTRUCT.motifnum);
bostypelist = unique(DATSTRUCT.bostype);

AllUnits_FRmat = {};
AllUnits_FRmat_t = {};
AllUnits_SylOnsets = {};
AllUnits_SylOffsets = {};
AllUnits_bostype = {};
AllUnits_motifnum = [];
AllUnits_exptnum = [];
AllUnits_unitnum = [];
AllUnits_bregions = {};

for i=1:length(indsgrp_unique)
    
    for mm=motifslist'
        for bb=1:length(bostypelist)
            bosthis = bostypelist(bb);
            
            indsthis = indsgrp==indsgrp_unique(i) & DATSTRUCT.motifnum==mm & DATSTRUCT.bostype==bosthis;
            if isempty(indsthis)
                continue
            end
            
            rendnum = DATSTRUCT.songrendnum(indsthis);
            spks = DATSTRUCT.Spktimes(indsthis);
            sylonsets = cell2mat(DATSTRUCT.Onstimes(indsthis));
            syloffsets = cell2mat(DATSTRUCT.Offtimes(indsthis));
            
            minmaxtime = DATSTRUCT.minmax_time(indsthis,:);
            
            % === for pplotting.
            birdname = unique(DATSTRUCT.birdname(indsthis));
            exptnum = unique(DATSTRUCT.exptnum(indsthis));
            unitnum = unique(DATSTRUCT.unitnum(indsthis));
            bregion = SummaryBOS.expt(exptnum).bregions{unitnum};
            motif = unique(DATSTRUCT.motifregexp(indsthis));
            bostype = PARAMS.BOSnames{strcmp(PARAMS.BOSbirdname, birdname)}{bosthis};
            
            %% ========== PLOT FOR THIS UNIT, RESPONSES OVER ALL RENDITIONS OF
            % THIS MOTIF
            assert(length(unique(minmaxtime(:,1)))==1, 'assumes all predur from same alignment');
            mintime = minmaxtime(1,1);
            maxtime = min(minmaxtime(:,2));
            
            
            [FRmat, t] = lt_neural_QUICK_Spk2FRmat(spks, mintime, maxtime);
            
            frmean = mean(FRmat,2);
            frsem = lt_sem(FRmat');
            
            
            % --- SYLS ONSET OFFSET
            on = median(sylonsets,1);
            off = median(syloffsets,1);
            
            
            % =========================== COLLECT
            AllUnits_FRmat = [AllUnits_FRmat; FRmat];
            AllUnits_FRmat_t = [AllUnits_FRmat_t; t];
            AllUnits_SylOnsets = [AllUnits_SylOnsets; sylonsets];
            AllUnits_SylOffsets = [AllUnits_SylOffsets; syloffsets];
            
            AllUnits_motifnum = [AllUnits_motifnum; mm];
            AllUnits_bostype = [AllUnits_bostype; bostype];
            AllUnits_exptnum = [AllUnits_exptnum; exptnum];
            AllUnits_unitnum = [AllUnits_unitnum; unitnum];
            
            AllUnits_bregions = [AllUnits_bregions; bregion];

        end
    end
end


%% ==================== [ANALYSIS - PLOT], MEAN RESPONSE OVER UNITS.
% ===== will automatically go thru all motifs.

exptlist = 1:3; % takes units from these experiments
exptlist = 4:7; % takes units from these experiments
bostype = 'fwd';
bregionthis = 'RA';

% --- make sure all are from same bird
birdname = unique({SummaryBOS.expt(exptlist).birdname});
assert(length(birdname)==1, 'all ext should be from one bird');


figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% #########################################
indsthis = strcmp(AllUnits_bregions, bregionthis) & ismember(AllUnits_exptnum, exptlist) ...
    & strcmp(AllUnits_bostype, bostype);

% ---- for each unit, get mean response to each mtoif
frmat = AllUnits_FRmat(indsthis);
frmat_t = AllUnits_FRmat_t(indsthis);
motifnum = AllUnits_motifnum(indsthis);
units = AllUnits_unitnum(indsthis);
expt = AllUnits_exptnum(indsthis);
units_unique = lt_tools_grp2idx({expt, units});

% --- GET MEAN RESPONSE FOR EACH UNIT/MOTIF
frmean = cellfun(@(x)mean(x,2), frmat, 'UniformOutput', 0);
frsem = cellfun(@(x)lt_sem(x'), frmat,  'UniformOutput', 0);
tmean = cellfun(@(x)mean(x,2), frmat_t, 'UniformOutput', 0);

% --- TAKE MEAN ACROSS UNITS
nummotifs = max(motifnum);
hsplots = [];
for mm=1:nummotifs
    frthis = frmean(motifnum==mm);
    tthis = tmean(motifnum==mm);
    
    ymean = mean(cell2mat(cellfun(@transpose, frthis, 'UniformOutput', 0)));
    ysem = lt_sem(cell2mat(cellfun(@transpose, frthis, 'UniformOutput', 0)));
    t = mean(cell2mat(cellfun(@transpose, tthis, 'UniformOutput', 0)));
    
    % ======== PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    motifthis = PARAMS.MotifsToGet{find(strcmp(PARAMS.MotifsToGet, birdname{1}))+1}{mm};
    title([birdname{1} '-expt' num2str(exptlist) '-' motifthis '[' bregionthis ']']);
   
    shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
    
    axis tight;
    lt_plot_zeroline;
end

linkaxes(hsplot, 'xy');

%% ==================== [ANALYSIS - PLOT], MEAN RESPONSE OVER UNITS.
% ===== will automatically go thru all motifs.
% NOTE: came back to this code, not sure what this sections id dpoig that
% is unique over previuos -- seems redundant?

% exptlist = 1:3; % takes units from these experiments
% bostype = 'fwd';
% bregionthis = 'LMAN';
% 
% % --- make sure all are from same bird
% assert(length(unique({SummaryBOS.expt(exptlist).birdname}))==1, 'all ext should be from one bird');
% 
% % ======== 
% indsthis = ismember(DATSTRUCT.exptnum, exptlist) & strcmp(DATSTRUCT.bostype_str, bostype) ...
%     & strcmp(DATSTRUCT.brainregion, bregionthis);
% 


%% ========================= [GOOD] - FIlter to get only units of interest.
% ============= input desired

exptlist_toget = 4:7; % takes units from these experiments
bostype_toget = 'fwd';
bregionthis_toget = 'RA';
% motif_toget = '(a)abhh';
motif_toget = 'jj(b)hh';


AllUnits = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS, PARAMS, ...
    exptlist_toget, bostype_toget, bregionthis_toget, motif_toget);


% ============ EXTRACT STATS ACROSS THESE UNITS 
[Stats] = lt_neural_BOS_StatsUnits(AllUnits, SummaryBOS, PARAMS);


% ======== PLOT
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

birdname = unique({SummaryBOS.expt(exptlist_toget).birdname}); birdname = birdname{1};
motifnum = unique(AllUnits.motifnum);
pcol = 'k';

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% hsplots = [hsplots; hsplot];
motifthis = PARAMS.MotifsToGet{find(strcmp(PARAMS.MotifsToGet, birdname))+1}{motifnum};
title([birdname '-expt' num2str(exptlist_toget) '-' motif_toget '[' bregionthis_toget ']']);

lt_neural_QUICK_PlotSylPatches(Stats.sylon_median, Stats.syloff_median, [0 max(Stats.frate_mean)], 0, pcol);
shadedErrorBar(Stats.frate_t, Stats.frate_mean, Stats.frate_sem, {'Color', pcol}, 1);

axis tight;
lt_plot_zeroline;
linkaxes(hsplot, 'xy');


%% ===================== [PLOT - GOOD] - each motif, mean across units
ExptlistList = {[1:3], [4:7]};
BosTypeList = {'fwd'};
BregionList = {'LMAN', 'RA'};
% MotifList = {'(a)abhh', '(j)jbhh'};
MotifList = {'aa(b)hh', 'jj(b)hh'};

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ============ iterate over: exptlist, bregion, bostype, motif
hsplots = [];
for j=1:length(ExptlistList)
    exptlist = ExptlistList{j};
    for jj=1:length(BosTypeList)
        bostype = BosTypeList{jj};
        for jjj=1:length(BregionList)
            bregion = BregionList{jjj};
            for k=1:length(MotifList)
                motif = MotifList{k};
                
                [AllUnits, birdname] = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS, PARAMS, ...
                    exptlist, bostype, bregion, motif);
                
                % ============ EXTRACT STATS ACROSS THESE UNITS
                [Stats] = lt_neural_BOS_StatsUnits(AllUnits, SummaryBOS, PARAMS);
                
                motifnum = unique(AllUnits.motifnum);
                pcol = 'k';
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots; hsplot];
                title([birdname '-expt' num2str(exptlist) '-' motif '[' bregion ']']);
                
                lt_neural_QUICK_PlotSylPatches(Stats.sylon_median, ...
                    Stats.syloff_median, [0 max(Stats.frate_mean)], 0, pcol);
                
                shadedErrorBar(Stats.frate_t, Stats.frate_mean, ...
                    Stats.frate_sem, {'Color', pcol}, 1);
                
                axis tight;
                lt_plot_zeroline;
            end
        end
    end
end
                linkaxes(hsplots, 'xy');



%% ====================== [ANALYSIS] - MOTIF MINUS MOTIF
%  MOTIF 1 MINUS MOTIF 2
motif1 = 'aa(b)hh';
motif2 = 'jj(b)hh';

ExptlistList = {[1:3], [4:7]};
BosTypeList = {'fwd'};
BregionList = {'LMAN', 'RA'};
% MotifList = {'(a)abhh', '(j)jbhh'};
MotifList = {'aa(b)hh', 'jj(b)hh'};

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ============ iterate over: exptlist, bregion, bostype, motif
hsplots = [];
for j=1:length(ExptlistList)
    exptlist = ExptlistList{j};
    for jj=1:length(BosTypeList)
        bostype = BosTypeList{jj};
        for jjj=1:length(BregionList)
            bregion = BregionList{jjj};
            
            pcol = 'k';
            [AllUnitsMotifDiff, birdname] = lt_neural_BOS_MotifDiff(DATSTRUCT, SummaryBOS, PARAMS, ...
                exptlist, bostype, bregion, motif1, motif2);
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots; hsplot];
            title([birdname '-expt' num2str(exptlist) '-' motif1 '-' motif2 '[' bregion ']']);
            
            % ======= plot difference for each unit
            y = cell2mat(cellfun(@transpose, AllUnitsMotifDiff.frdiff, 'UniformOutput', 0));
            t = AllUnitsMotifDiff.t;
            plot(t{1}, y', 'Color', [0.7 0.7 0.7]);
            
            % ======= mean difference across units
            ymean = mean(y,1);
            ysem = lt_sem(y);
            shadedErrorBar(t{1}, ymean, ysem, {'Color', pcol},1);
            
            % ====== plot syl patches
            [AllUnits] = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS, PARAMS, ...
                exptlist, bostype, bregion, motif1);
            [Stats] = lt_neural_BOS_StatsUnits(AllUnits, SummaryBOS, PARAMS);
            lt_neural_QUICK_PlotSylPatches(Stats.sylon_median, ...
                Stats.syloff_median, [0 max(Stats.frate_mean)], 0, 'r');

            [AllUnits] = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS, PARAMS, ...
                exptlist, bostype, bregion, motif2);
            [Stats] = lt_neural_BOS_StatsUnits(AllUnits, SummaryBOS, PARAMS);
            lt_neural_QUICK_PlotSylPatches(Stats.sylon_median, ...
                Stats.syloff_median, [0 max(Stats.frate_mean)], 0, 'b');

            axis tight;
            lt_plot_zeroline;
        end
    end
end
linkaxes(hsplots, 'xy');



