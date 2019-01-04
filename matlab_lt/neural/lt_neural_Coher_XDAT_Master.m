%%  =========== 2) SUMMARY PLOT ACROSS MULITPLE MOTIFPLOTS

disp('NOTE: this assumes that each saved dataset saved after taking across channel means...');
disp('ALSO assumes that all have same time and freq bins, and those bins are in current PARAMS in worksp[ace.');
pause

savedir = '/bluejay0/bluejay2/lucas/analyses/neural/COHERENCE/SaveMotifDat';
MotifPlotList = {'15Dec2018_1321', '20Dec2018_0951', '03Jan2019_0153', '03Jan2019_2009'};
MotifPlotDatType = {'data', 'data', 'negcontrol', 'data'};

% ================= FILTER WHAT CASES YOU WANT TO EXTRACT
% --- to get specific switch types. ... [is done in addition to above
% fitlers]
% swtoget = {}; % passes if matches ANY of these
swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
% for a given switch did not have a previous switch on the same day
firstswitch=1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitch);
% indtoget_b_e_s = [];

% ========================== ANY THINGS TO SKIP?
if ~isempty(indtoget_b_e_s)
%     indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1 2 3]), [3 1 7], 'rows'), :) = []; % this expt no learning.
%     indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1 2 3]), [3 1 11], 'rows'), :) = []; % this expt need to remove off target syl first...
    % indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1]), [3], 'rows'), :) = []; % this bird did not have RAoutside recordings

    indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1]), [2], 'rows'), :) = []; % this bird did not have RAoutside recordings
end


[XDAT] = lt_neural_Coher_XDATSET_Extract(savedir, indtoget_b_e_s, SwitchStruct);


% ====== break out XDAT (i.e. fields into the variables)

tmp = fieldnames(XDAT);
for i=1:length(tmp)
    eval([tmp{i} ' = XDAT.(tmp{i});']);
end


%% ============ extract dat type

all_dattype = cellfun(@(x)(MotifPlotDatType{strcmp(MotifPlotList, x)}), all_savemarker, 'UniformOutput', 0);



%% =========== REDEFINE T AND F WINDOW FOR COHERENCE SCALAR?
tbins = [-0.07 -0.03];
% fbins = [20 35];
fbins = [20 35];

% ============== RUNS
inds_t = PARAMS.tbins>tbins(1) & PARAMS.tbins<tbins(2);
inds_f = PARAMS.ffbins>fbins(1) & PARAMS.ffbins<fbins(2);

for i=1:length(all_data)
   
tmp = squeeze(mean(mean(all_data(i).cohmat_diff(inds_t, inds_f, :),1),2));
all_data(i).cohscal_diff = tmp; 
    
end

%% =========== 2b) [SUMMARY PLOT] - ALONG MOTIF
close all;
subtractAcrossSylMean = 1;

[indsgrp, indsgrp_unique] = lt_tools_grp2idx({all_bname, all_ename, all_swnum});

figcount=1;
subplotrows=5;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ======= 1) FOR EACH SWITCH, OVERLAY DATA AND CONTROL
for i=1:length(indsgrp_unique)
    indsthis = find(indsgrp == indsgrp_unique(i));
    
    bname = all_bname{indsthis(1)};
    ename = all_ename{indsthis(1)};
    swnum = all_swnum(indsthis(1));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel('cohscal (WN - base) (rd = WN)');
    for ii=1:length(indsthis)
        
        datthis = all_data(indsthis(ii));
        savemark = all_savemarker{indsthis(ii)};
        dattype = MotifPlotDatType{strcmp(MotifPlotList, savemark)};
        
        if strcmp(dattype, 'data')
            pcol = 'r';
        elseif strcmp(dattype, 'negcontrol')
            pcol = [0.6 0.6 0.6];
        else
            asdfasfd;
        end
        
        % --- plot coh scalar (WN minsu base)
        if subtractAcrossSylMean==1
        y = datthis.cohscal_diff - mean(datthis.cohscal_diff);
        
            ylabel('cohscal (WN - base) (rd = WN) [subtract across syl mean]');

        else
            y = datthis.cohscal_diff;
    ylabel('cohscal (WN - base) (rd = WN)');
        end
        
        plot(datthis.motifID, y, '-o', 'Color', pcol);
        
        YLIM = ylim;
        
        % --- plot whehter is targ/same
        plot(datthis.motifID(datthis.istarg==1), YLIM(1)+0.02, 'r^');
        if any(datthis.issame)
            plot(datthis.motifID(datthis.issame==1), YLIM(1)+0.02, 'b^');
        end
    end
    
    % --- figure out what the motifs names are
    [~, motiflist_out, ~] = ...
        lt_neural_QUICK_MotifID(bname, '');
    set(gca, 'XTick', 1:length(motiflist_out), 'XTickLabel', motiflist_out);
    rotateXLabels(gca, 90);
    
    xlim([0 length(motiflist_out)+1]);
    lt_plot_zeroline;
end



%% =========== 2b) [SUMMARY PLOT] - COHEROGRAMS [TARG, SAME, DIFF, AND TARG-REST]
close all;

[indsgrp, indsgrp_unique] = lt_tools_grp2idx({all_bname, all_ename, all_swnum});

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ======= 1) FOR EACH SWITCH, OVERLAY DATA AND CONTROL
for i=1:length(indsgrp_unique)
    indsthis = find(indsgrp == indsgrp_unique(i));
    
    bname = all_bname{indsthis(1)};
    ename = all_ename{indsthis(1)};
    swnum = all_swnum(indsthis(1));
    
    for ii=1:length(indsthis) % diff save structures
        
        datthis = all_data(indsthis(ii));
        savemark = all_savemarker{indsthis(ii)};
        dattype = MotifPlotDatType{strcmp(MotifPlotList, savemark)};
        %
        %         if strcmp(dattype, 'data')
        %             pcol = 'r';
        %         elseif strcmp(dattype, 'negcontrol')
        %             pcol = [0.6 0.6 0.6];
        %         else
        %             asdfasfd;
        %         end
        
        % === plot mean coherence matrix for target, same, and diff
        cohmat_targ = mean(datthis.cohmat_diff(:,:, datthis.istarg==1),3);
        cohmat_same = mean(datthis.cohmat_diff(:,:, datthis.istarg==0 & datthis.issame==1),3);
        cohmat_diff = mean(datthis.cohmat_diff(:,:, datthis.istarg==0 & datthis.issame==0),3);
        cohmat_targminrest = cohmat_targ - ...
            mean(datthis.cohmat_diff(:,:, datthis.istarg==0),3);
        
        % ================ TARG
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title({[bname '-' ename], ['-sw' num2str(swnum) ', ' dattype]});
        ylabel('TARG');
        %         ylabel('cohscal (WN - base) (rd = WN)');
        lt_neural_Coher_Plot(cohmat_targ, PARAMS.tbins, PARAMS.ffbins, 1, '', ...
            [-0.15 0.15]);
        colorbar('East');
        
        
        % ================ SAME
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title({[bname '-' ename], ['-sw' num2str(swnum) ', ' dattype]});
        ylabel('SAME');
        %         ylabel('cohscal (WN - base) (rd = WN)');
        lt_neural_Coher_Plot(cohmat_same, PARAMS.tbins, PARAMS.ffbins, 1, '', ...
            [-0.15 0.15]);
        colorbar('East');
        
        
        % ================ DIFF
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title({[bname '-' ename], ['-sw' num2str(swnum) ', ' dattype]});
        ylabel('DIFF');
        %         ylabel('cohscal (WN - base) (rd = WN)');
        lt_neural_Coher_Plot(cohmat_diff, PARAMS.tbins, PARAMS.ffbins, 1, '', ...
            [-0.15 0.15]);
        colorbar('East');
        
        
        % ================= TARG MINUS REST
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title({[bname '-' ename], ['-sw' num2str(swnum) ', ' dattype]});
        ylabel('TARG - REST');
        %         ylabel('cohscal (WN - base) (rd = WN)');
        lt_neural_Coher_Plot(cohmat_targminrest, PARAMS.tbins, PARAMS.ffbins, 1, '', ...
            [-0.15 0.15]);
        colorbar('East');
        
        
    end
    
end


%% =========== 2b) [SUMMARY PLOT] - COHEROGRAMS [COMPARING DAT VS. CONTROL]
close all;

[indsgrp, indsgrp_unique] = lt_tools_grp2idx({all_bname, all_ename, all_swnum});

figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ======= 1) FOR EACH SWITCH, OVERLAY DATA AND CONTROL
Yallmat = {}; % -- collect all summary mats
for i=1:length(indsgrp_unique)
    indsthis = find(indsgrp == indsgrp_unique(i));

    bname = all_bname{indsthis(1)};
    ename = all_ename{indsthis(1)};
    swnum = all_swnum(indsthis(1));
            
    % ==== ONLY CONTINUE IF THIS EXPERIMENT HAS BOTH DAT AND CONTROLS
    if length(unique(all_dattype(indsthis)))<2
        % then doesn't have both dat and control...
        continue
    else
        assert(length(unique(all_dattype(indsthis)))==2, 'multiple dataset of same thing...');
    end
    
    Y = cell(1, length(indsthis));
    for ii=1:length(indsthis) % diff save structures

        datthis = all_data(indsthis(ii));
        savemark = all_savemarker{indsthis(ii)};
        dattype = MotifPlotDatType{strcmp(MotifPlotList, savemark)};
        
        % === plot mean coherence matrix for target, same, and diff
        cohmat_targ = mean(datthis.cohmat_diff(:,:, datthis.istarg==1),3);
        cohmat_same = mean(datthis.cohmat_diff(:,:, datthis.istarg==0 & datthis.issame==1),3);
        cohmat_diff = mean(datthis.cohmat_diff(:,:, datthis.istarg==0 & datthis.issame==0),3);
        cohmat_targminrest = cohmat_targ - ...
            mean(datthis.cohmat_diff(:,:, datthis.istarg==0),3);
        
        if strcmp(dattype, 'data')
            Y{2} = cohmat_targminrest;
        elseif strcmp(dattype, 'negcontrol')
            Y{1} = cohmat_targminrest;
        end
    end
    Ythis = Y{2} - Y{1};
    Yallmat = [Yallmat; Ythis];

    % ================= 
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title({[bname '-' ename], ['-sw' num2str(swnum) ', ' dattype]});
        ylabel('DAT - Control (of TARG - REST)');
        %         ylabel('cohscal (WN - base) (rd = WN)');
        lt_neural_Coher_Plot(Ythis, PARAMS.tbins, PARAMS.ffbins, 1, '', ...
            [-0.15 0.15]);
        colorbar('East');
        ylim([20 100]);
end

% === summary (mean of all)
Yallmat = lt_neural_Coher_Cell2Mat(Yallmat); 
% mean(Yallmat,3);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('mean across expt');
        ylabel('DAT - Control (of TARG - REST)');
        lt_neural_Coher_Plot(Yallmat, PARAMS.tbins, PARAMS.ffbins, 1, '', ...
            [-0.15 0.15]);
        colorbar('East');
        ylim([20 100]);


%% ===== [SCALAR] - EACH EXPERIMENT ONE DOT
% 1) DAT MINUS CONTROL
% 2) DAT (TARG MINUS REST)

close all;

% % =====================================
% % --- to get specific switch types. ... [is done in addition to above
% % fitlers]
% % swtoget = {}; % passes if matches ANY of these
% swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% % swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
% firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
% % for a given switch did not have a previous switch on the same day
% indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
%     firstswitchfortarget_withinday);
% % indtoget_b_e_s = [];
% 
% % ========================== ANY THINGS TO SKIP?
% if (1)
%     indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1 2 3]), [3 1 7], 'rows'), :) = []; % this expt no learning.
%     indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1 2 3]), [3 1 11], 'rows'), :) = []; % this expt need to remove off target syl first...
%     indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1]), [2], 'rows'), :) = []; % this bird did not have RAoutside recordings
%     % indtoget_b_e_s(ismember(indtoget_b_e_s(:,[1]), [3], 'rows'), :) = []; % this bird did not have RAoutside recordings
% end
% 
% 
% ===================================== 1) EXTRACT ONE SCALAR EACH GROUP
[indsgrp, indsgrp_unique] = lt_tools_grp2idx({all_bname, all_ename, all_swnum});

% ========= 1) COLLECT ALL DATA (each swtich get individual scalars)
allgrp_grpind = [];
allgrp_cohtargminusrest = [];
allgrp_dattype = []; % 0=control; 1=data
% --- bird/e/s correspond to what is in switchstruct.
allgrp_bnum = [];
allgrp_enum = [];
allgrp_swnum = [];

for i=1:length(indsgrp_unique)
    indsthis = find(indsgrp == indsgrp_unique(i));
    
    bname = all_bname{indsthis(1)};
    ename = all_ename{indsthis(1)};
    swnum = all_swnum(indsthis(1));
    
    %     if strcmp(bname, 'wh44wh39') & strcmp(ename, 'RALMANlearn1') & swnum==2
    %         keyboard
    %     end
    
    bnum = find(strcmp({SwitchStruct.bird.birdname}, bname));
    enum = find(strcmp({SwitchStruct.bird(bnum).exptnum.exptname}, ename));
    
%     % ==== decide whether to skip
%     if ~isempty(indtoget_b_e_s)
%         if ~ismember([bnum enum swnum], indtoget_b_e_s, 'rows')
%             continue
%         end
%     end
    
    for ii=1:length(indsthis)
        
        datthis = all_data(indsthis(ii));
        savemark = all_savemarker{indsthis(ii)};
        dattype = MotifPlotDatType{strcmp(MotifPlotList, savemark)};
        
        % ===== EXTRACT TARG MINUS REST
        cohtarg = mean(datthis.cohscal_diff(datthis.istarg==1));
        cohsame = mean(datthis.cohscal_diff(datthis.istarg==0 & datthis.issame==1));
        cohdiff = mean(datthis.cohscal_diff(datthis.istarg==0 & datthis.issame==0));
        cohtargminusrest = cohtarg - mean(datthis.cohscal_diff(datthis.istarg==0));
%         cohtargminusrest = cohtarg - cohdiff;
        
        allgrp_grpind = [allgrp_grpind; indsgrp_unique(i)];
        allgrp_cohtargminusrest = [allgrp_cohtargminusrest; cohtargminusrest];
        if strcmp(dattype, 'data')
            allgrp_dattype = [allgrp_dattype; 1];
        elseif strcmp(dattype, 'negcontrol')
            allgrp_dattype = [allgrp_dattype; 0];
        end
        
        
        % ============== what bird/expt?
        allgrp_bnum = [allgrp_bnum; bnum];
        allgrp_enum = [allgrp_enum; enum];
        allgrp_swnum = [allgrp_swnum; swnum];
    end
end


% ============= PLOT [EACH EXPERIMENT, DATA VS. CONTROL]
lt_figure; hold on;
xlabel('DATA --- CONTROL');
ylabel('coh (targ minus rest');
YY = [];
plotcols = lt_make_plot_colors(max(allgrp_bnum), 0, 0);

for i=unique(allgrp_grpind)'
    
    Y = nan(1,2);
    pcol = plotcols{unique(allgrp_bnum(allgrp_grpind==i))};
    
    % -- data
    indsthis = allgrp_grpind==i & allgrp_dattype==1;
    coh = mean(allgrp_cohtargminusrest(indsthis));
    Y(1) = coh;
    
    
    % -- control
    indsthis = allgrp_grpind==i & allgrp_dattype==0;
    coh = mean(allgrp_cohtargminusrest(indsthis));
    Y(2) = coh;
    
    % ---- PLOT
    x = [1 2];
    plot(x, Y, '-o', 'Color', pcol);
    
    bname = SwitchStruct.bird(unique(allgrp_bnum(allgrp_grpind==i))).birdname;
    ename = SwitchStruct.bird(unique(allgrp_bnum(allgrp_grpind==i))).exptnum(unique(allgrp_enum(allgrp_grpind==i))).exptname;
    swnum = unique(allgrp_swnum(allgrp_grpind==i));
    lt_plot_text(x(end)+0.2, Y(end), [bname '-' ename '-sw' num2str(swnum)], 'm', 8);
    
    YY = [YY; Y];
end
xlim([0 3]);
lt_plot_zeroline;

lt_plot([1 2]+0.2, nanmean(YY,1), {'Errors', lt_sem(YY)});

[~, p] = ttest(YY(:,1), YY(:,2));
p = signrank(YY(:,1), YY(:,2));
lt_plot_pvalue(p, 'srank',1);


%%  ##############################################################
%% ========= PLOT SUMMARY FIGURES (ONLY DATA, NOT CONTROLS)
close all;

[indsgrp, indsgrp_unique] = lt_tools_grp2idx({all_bname, all_ename, all_swnum});

% ========= 1) COLLECT ALL DATA (each swtich get individual scalars)
allgrp_grpind = [];
allgrp_Y = []; % targ, same, diff, rest
% allgrp_cohtarg = [];
% allgrp_dattype = []; % 0=control; 1=data
% --- bird/e/s correspond to what is in switchstruct.
allgrp_bnum = [];
allgrp_enum = [];
allgrp_swnum = [];

for i=1:length(indsgrp_unique)
    indsthis = find(indsgrp == indsgrp_unique(i) & strcmp(all_dattype, 'data'));
    
    if isempty(indsthis)
        disp('x');
        continue
    end
    bname = all_bname{indsthis(1)};
    ename = all_ename{indsthis(1)};
    swnum = all_swnum(indsthis(1));
    
    %     if strcmp(bname, 'wh44wh39') & strcmp(ename, 'RALMANlearn1') & swnum==2
    %         keyboard
    %     end
    
    bnum = find(strcmp({SwitchStruct.bird.birdname}, bname));
    enum = find(strcmp({SwitchStruct.bird(bnum).exptnum.exptname}, ename));
    
%     % ==== decide whether to skip
%     if ~isempty(indtoget_b_e_s)
%         if ~ismember([bnum enum swnum], indtoget_b_e_s, 'rows')
%             continue
%         end
%     end
    
    for ii=1:length(indsthis)
        
        datthis = all_data(indsthis(ii));
        savemark = all_savemarker{indsthis(ii)};
        dattype = MotifPlotDatType{strcmp(MotifPlotList, savemark)};
        assert(strcmp(dattype, 'data'));
        
        % ===== EXTRACT TARG MINUS REST
        cohtarg = mean(datthis.cohscal_diff(datthis.istarg==1));
        cohsame = mean(datthis.cohscal_diff(datthis.istarg==0 & datthis.issame==1));
        cohdiff = mean(datthis.cohscal_diff(datthis.istarg==0 & datthis.issame==0));
        cohrest = mean(datthis.cohscal_diff(datthis.istarg==0));
%         cohtargminusrest = cohtarg - cohdiff;
        
        allgrp_grpind = [allgrp_grpind; indsgrp_unique(i)];
        allgrp_Y = [allgrp_Y; [cohtarg cohsame cohdiff cohrest]];
        
        if strcmp(dattype, 'data')
            allgrp_dattype = [allgrp_dattype; 1];
        elseif strcmp(dattype, 'negcontrol')
            allgrp_dattype = [allgrp_dattype; 0];
        end
        
        
        % ============== what bird/expt?
        allgrp_bnum = [allgrp_bnum; bnum];
        allgrp_enum = [allgrp_enum; enum];
        allgrp_swnum = [allgrp_swnum; swnum];
    end
end


% =========== PLOTS
lt_figure; hold on;

% === 1) 
lt_subplot(3,2,1); hold on;
xlabel('targ -- rest');
ylabel('coh (wn - base)');

colsthis = [1 4];
indsthis = all(~isnan(allgrp_Y(:, colsthis))');
y = allgrp_Y(indsthis,colsthis);
x = colsthis;
plot(x, y, '-ok');
xlim([0 5]);
if length(colsthis)==2
   p = signrank(y(:,1), y(:,2));
   lt_plot_pvalue(p, 'srank', 1);
elseif length(colsthis)==3
   p = signrank(y(:,1), y(:,2));
   lt_plot_pvalue(p, 'srank(1vs2)', 1);
   p = signrank(y(:,1), y(:,3));
   lt_plot_pvalue(p, 'srank(1vs3)', 2);
end

% === 2) 
lt_subplot(3,2,2); hold on;
xlabel('targ -- same - diff');
ylabel('coh (wn - base)');

colsthis = [1 2 3];
indsthis = all(~isnan(allgrp_Y(:, colsthis))');
y = allgrp_Y(indsthis,colsthis);
x = colsthis;
plot(x, y, '-ok');
xlim([0 5]);
if length(colsthis)==2
   p = signrank(y(:,1), y(:,2));
   lt_plot_pvalue(p, 'srank', 1);
elseif length(colsthis)==3
   p = signrank(y(:,1), y(:,2));
   lt_plot_pvalue(p, 'srank(1vs2)', 1);
   p = signrank(y(:,1), y(:,3));
   lt_plot_pvalue(p, 'srank(1vs3)', 2);
end


