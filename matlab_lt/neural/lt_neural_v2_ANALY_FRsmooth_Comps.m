function OUTDAT = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, shuffSylType, ...
    epochtoplot, plotOn, shuffonlynontargs)
%% also collects params for learning.
%% lt 9/12/18 - comparisons between syl types
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
% shuffSylType = 0; % i.e. do shuffle of syl types within each neuron?
% plotOn =1;

if ~exist('shuffonlynontargs', 'var')
    shuffonlynontargs = 0;
end
%% DO SHUFFLE?
% ===============
maxneur = max(OUTDAT.All_neurnum);
numbirds = length(SwitchStruct.bird);

% ======================== SHUFFLE SYL TYPE? % within each neuron
if shuffSylType==1
    for i=1:numbirds
        numexpts = length(SwitchStruct.bird(i).exptnum);
        for ii=1:numexpts
            numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
            for ss = 1:numswitch
                % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
                for nn=1:maxneur
                    
                    if shuffonlynontargs==1
                        % -- then only get nontargs
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                        & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==0);
                    elseif shuffonlynontargs==0
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                        & OUTDAT.All_neurnum==nn);
                    end
                    
                    if ~any(indsthis)
                        continue
                    end
                    
                    % ================ SHUFFLE INDS AND REPLACE SYL TYPES
                    indstmp = indsthis(randperm(length(indsthis)));
%                     disp([num2str(indstmp') '-' num2str(indsthis')]);
                    OUTDAT.All_issame(indsthis) = OUTDAT.All_issame(indstmp);
                    OUTDAT.All_istarg(indsthis) = OUTDAT.All_istarg(indstmp);
                    OUTDAT.All_motifnum(indsthis) = OUTDAT.All_motifnum(indstmp);
%                     OUTDAT.All_FF(indsthis, :) = OUTDAT.All_FF(indstmp, :);
%                     OUTDAT.All_FF_t(indsthis, :) = OUTDAT.All_FF_t(indstmp, :);
                    
                end
            end
        end
    end
end


%% ====================== collect data
AllOnlyMinusBase_FRsmooth = cell(size(OUTDAT.All_FRsmooth,1),1); % subtract baseline
% === absoulte balues
AllMinusAll_FRsmooth = cell(size(OUTDAT.All_FRsmooth,1),1); % subtract baseline, then abs(subtract global mean)
AllMinusAllMinusDiff_FRsmooth = cell(size(OUTDAT.All_FRsmooth,1),1); % above, one step further, subtract mean of diff types.
AllOnlyMinusDiff_FRsmooth = cell(size(OUTDAT.All_FRsmooth,1),1); % subtract base, then abs dev from mean of diff type
% --- not absolute balues
AllDevDiff_NotAbs = cell(size(OUTDAT.All_FRsmooth,1),1); % subtract base, then dev from mean of diff typ


for i=1:numbirds
    
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if ~any(indsthis)
                    continue
                end
                
                % ================ COLLECT FRATE
                Sameall = OUTDAT.All_issame(indsthis);
                Targall = OUTDAT.All_istarg(indsthis);
                Yall = [];
                %                 tbin = OUTDAT.All_FRsmooth_t(find(indsthis,1, 'first'));
                for j=indsthis'
                    %
                    %                     [FRmeanAll, FRsemAll, tbin] = fn_subtractbase(OUTDAT, j, prctile_divs, usepercent);
                    %
                    Yall = [Yall OUTDAT.AllMinusBase_FRmeanAll{j}{epochtoplot}];
                end
                
                % ===== to save only minus base
                Yall_orig = Yall;
                
                % ############################# DEVIATION FROM DIFF TYPE
                YdevFromDiff = Yall - mean(Yall(:, Sameall==0 & Targall==0), 2);
                
                % ############################# ABSOLUTE VALUES
                % === for each one, get absolute value deviation from mean
                % of different types
                Yall_onlyminusdiff = abs(Yall - mean(Yall(:, Sameall==0 & Targall==0), 2));
                
                % ==== for each one, get absolute value deviation from
                % global mean
                ymean_all = mean(Yall,2);
                Yall = abs(Yall - ymean_all);
                
                
                % ==== subtract mean of diff type
                ymean_difftype = mean(Yall(:, Sameall==0 & Targall==0), 2);
                Yall_minusDiff = Yall - ymean_difftype;
                
                % ==== STORE DATA
                for j=1:size(Yall,2)
                    
                    % ===== not absoute values
                    AllOnlyMinusBase_FRsmooth{indsthis(j)} = Yall_orig(:,j)';
                    AllDevDiff_NotAbs{indsthis(j)} = YdevFromDiff(:,j);
                    
                    % === absolute values
                    AllMinusAll_FRsmooth{indsthis(j)} = Yall(:,j)';
                    AllMinusAllMinusDiff_FRsmooth{indsthis(j)} = Yall_minusDiff(:,j)';
                    AllOnlyMinusDiff_FRsmooth{indsthis(j)} = Yall_onlyminusdiff(:,j)';
                end
                
            end
        end
    end
end
OUTDAT.AllDevDiff_NotAbs = AllDevDiff_NotAbs;
OUTDAT.AllOnlyMinusBase_FRsmooth = AllOnlyMinusBase_FRsmooth;
OUTDAT.AllMinusAll_FRsmooth = AllMinusAll_FRsmooth;
OUTDAT.AllMinusAllMinusDiff_FRsmooth = AllMinusAllMinusDiff_FRsmooth;
OUTDAT.AllOnlyMinusDiff_FRsmooth = AllOnlyMinusDiff_FRsmooth;


%%
if plotOn==1
    % ####################################### PLOT SUMMARY
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ============ 1) DEVIATION FROM BASE
    % ----------- TARG
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('targ, all');
    ylabel('deviation from base)');
    pcol = 'k';
    indsthis = OUTDAT.All_istarg==1;
    
    ymat = cell2mat(OUTDAT.AllOnlyMinusBase_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    % ----------- SAME
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same, all');
    ylabel('deviation from base)');
    pcol = 'b';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==1;
    
    ymat = cell2mat(OUTDAT.AllOnlyMinusBase_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    % ----------- DIFF
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff, all');
    ylabel('deviation from base)');
    pcol = 'r';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==0;
    
    ymat = cell2mat(OUTDAT.AllOnlyMinusBase_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    
    
    % ============ 1) SUBTRACT MEAN OF DIFF
    % ----------- TARG
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('targ, all');
    ylabel('abs frate minsu base, then minus mean(difftype)');
    pcol = 'k';
    indsthis = OUTDAT.All_istarg==1;
    
    ymat = cell2mat(OUTDAT.AllOnlyMinusDiff_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    % ----------- SAME
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same, all');
    ylabel('abs frate minsu base, then minus mean(difftype)');
    pcol = 'b';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==1;
    
    ymat = cell2mat(OUTDAT.AllOnlyMinusDiff_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    % ----------- DIFF
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff, all');
    ylabel('abs frate minsu base, then minus mean(difftype)');
    pcol = 'r';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==0;
    
    ymat = cell2mat(OUTDAT.AllOnlyMinusDiff_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    
    
    % ========== 1) all trials + global mean (targ)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('targ, all');
    ylabel('abs frate minsu base, then minus global');
    pcol = 'k';
    indsthis = OUTDAT.All_istarg==1;
    
    ymat = cell2mat(OUTDAT.AllMinusAll_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    % ========== 1) all trials + global mean (targ)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same, all');
    pcol = 'b';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==1;
    
    ymat = cell2mat(OUTDAT.AllMinusAll_FRsmooth(indsthis));
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    % ========== 1) all trials + global mean (targ)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff, all');
    pcol = 'r';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==0;
    
    ymat = cell2mat(OUTDAT.AllMinusAll_FRsmooth(indsthis));
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    
    
    
    % %%%%%%%%%%%%%%%%%% SUBTRACT DIFF TYPES
    % -------------- TARG
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('targ, all');
    ylabel('minus base, abs(minus mean all), minus mean diff');
    pcol = 'k';
    indsthis = OUTDAT.All_istarg==1;
    
    ymat = cell2mat(OUTDAT.AllMinusAllMinusDiff_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    lt_plot_zeroline;
    
    % ------------- SAME
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('same, all');
    ylabel('minus base, abs(minus mean all), minus mean diff');
    pcol = 'b';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==1;
    
    ymat = cell2mat(OUTDAT.AllMinusAllMinusDiff_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    lt_plot_zeroline;
    
    % ------------- DIFF
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff, all');
    ylabel('minus base, abs(minus mean all), minus mean diff');
    pcol = 'r';
    indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==0;
    
    ymat = cell2mat(OUTDAT.AllMinusAllMinusDiff_FRsmooth(indsthis));
    
    t = OUTDAT.All_FRsmooth_t{1};
    plot(t, ymat', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(t, mean(ymat,1), lt_sem(ymat), {'Color', pcol}, 1);
    lt_plot_zeroline;
    
end


%% ========= PLOT DEVIATIONS FROM DIFF TYPE SEPARATELY DEPEND ON BASELINE CORR WITH PITCH



