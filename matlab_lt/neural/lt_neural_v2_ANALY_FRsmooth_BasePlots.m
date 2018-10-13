function lt_neural_v2_ANALY_FRsmooth_BasePlots(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, epochtoplot, plotDevFromBase, plotNormMeasures)
%% ####################### GETTING DEVIATION FROM BASELINE SMOOTHED FR

% prctile_divs = [33 66 100]; % percentiles to divide up data by
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
% prctile_divs = [50 100]; % percentiles to divide up data by
% epochtoplot = 2; % i.e. out of the epochs decided by prctile_divs

% usepercent = 0; % if 1, then gets fr percent deviation from baseline. otherwise uses hz diff
% nbasetime = 60; % 60 minutes bnefore first WN trial is min time for baseline
% nbasetime = []; % 60 minutes bnefore first WN trial is min time for baseline

numbirds = length(SwitchStruct.bird);
maxneur = max(OUTDAT.All_neurnum);

%% ############# PLOTS OF RAW DATA
%%
% nbase = 20; % get last 20 renditions
if plotDevFromBase==1


for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        
        % ######################## FIRST, PLOT LEARNING FOR ALL SYLS
        bname = MOTIFSTATS_Compiled.birds(i).birdname;
        ename = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        MeanSubtract = 1;
        BirdsToPlot = {bname};
        ExptToPlot = {ename}; % expt and bird must intersect.
        lt_neural_v2_ANALY_LearnAllSylPlot(MOTIFSTATS_Compiled, ...
            MeanSubtract, BirdsToPlot,ExptToPlot);
        
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            figcount=1;
            subplotrows=8;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            for nn=1:maxneur
                
                if ~any(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                        & OUTDAT.All_neurnum==nn)
                    continue
                end
                
                
                % ======= targ
                istarg_this = 1;
                issame_this = 1;
                pcol = 'k';
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['targ']);
                ylabel({[bname], [ename] ['sw' num2str(ss) ', neur' num2str(nn)]});
                ylabel('raw dev from baseline');
                %     ylabel(['neur ' num2str(nn) ', sw' num2str(ss)]);
                hsplots = [hsplots hsplot];
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==istarg_this & OUTDAT.All_issame==issame_this);
                
                yall = [];
                for j=indsthis'
                    
                    FRmeanAll = OUTDAT.AllMinusBase_FRmeanAll{j};
                    FRsemAll = OUTDAT.AllMinusBase_FRsemAll{j};
                    tbin = OUTDAT.AllMinusBase_tbinAll{j};
                    %                         [FRmeanAll, FRsemAll, tbin] = fn_subtractbase(OUTDAT, j, prctile_divs, usepercent, nbasetime);
                    
                    plot(tbin, FRmeanAll{epochtoplot}, '-', 'Color', pcol);
                    shadedErrorBar(tbin, FRmeanAll{epochtoplot}, FRsemAll{epochtoplot}, ...
                        {'Color', pcol},1);
                    
                    yall = [yall FRmeanAll{epochtoplot}];
                end
                plot(tbin, mean(yall,2), 'Color', pcol, 'LineWidth', 3);
                % ----- format stuff
                axis tight;
                line([0 0], ylim)
                lt_plot_zeroline;
                
                
                
                % ========= SAME
                istarg_this = 0;
                issame_this = 1;
                pcol = 'b';
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['same']);
                hsplots = [hsplots hsplot];
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==istarg_this & OUTDAT.All_issame==issame_this);
                
                if any(indsthis)
                    yall = [];
                    for j=indsthis'
                        
                        FRmeanAll = OUTDAT.AllMinusBase_FRmeanAll{j};
                        FRsemAll = OUTDAT.AllMinusBase_FRsemAll{j};
                        tbin = OUTDAT.AllMinusBase_tbinAll{j};
                        
                        plot(tbin, FRmeanAll{epochtoplot}, '-', 'Color', pcol);
                        shadedErrorBar(tbin, FRmeanAll{epochtoplot}, FRsemAll{epochtoplot}, ...
                            {'Color', pcol},1);
                        
                        yall = [yall FRmeanAll{epochtoplot}];
                    end
                    plot(tbin, mean(yall,2), 'Color', pcol, 'LineWidth', 3);
                    % ----- format stuff
                    axis tight;
                    line([0 0], ylim)
                    lt_plot_zeroline;
                end
                
                
                % ========= DIFF
                istarg_this = 0;
                issame_this = 0;
                pcol = 'r';
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['diff']);
                hsplots = [hsplots hsplot];
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==istarg_this & OUTDAT.All_issame==issame_this);
                
                if any(indsthis)
                    
                    yall = [];
                    for j=indsthis'
                        
                        
                        FRmeanAll = OUTDAT.AllMinusBase_FRmeanAll{j};
                        FRsemAll = OUTDAT.AllMinusBase_FRsemAll{j};
                        tbin = OUTDAT.AllMinusBase_tbinAll{j};
                        
                        plot(tbin, FRmeanAll{epochtoplot}, '-', 'Color', pcol);
                        shadedErrorBar(tbin, FRmeanAll{epochtoplot}, FRsemAll{epochtoplot}, ...
                            {'Color', pcol},1);
                        
                        yall = [yall FRmeanAll{epochtoplot}];
                    end
                    plot(tbin, mean(yall,2), 'Color', pcol, 'LineWidth', 3);
                    % ----- format stuff
                    axis tight;
                    line([0 0], ylim)
                    lt_plot_zeroline;
                end
                
                % ====================== format
                linkaxes(hsplots, 'xy');
                
            end
        end
    end
end
end

%% ############ ACCOUNT FOR DRIFT (take deviation from mean across al syls)
if plotNormMeasures==1
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        
        %         % ######################## FIRST, PLOT LEARNING FOR ALL SYLS
        bname = MOTIFSTATS_Compiled.birds(i).birdname;
        ename = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        %         MeanSubtract = 1;
        %         BirdsToPlot = {bname};
        %         ExptToPlot = {ename}; % expt and bird must intersect.
        %         lt_neural_v2_ANALY_LearnAllSylPlot(MOTIFSTATS_Compiled, ...
        %             MeanSubtract, BirdsToPlot,ExptToPlot);
        
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            figcount=1;
            subplotrows=8;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            for nn=1:maxneur
                
                if ~any(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                        & OUTDAT.All_neurnum==nn)
                    continue
                end
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                Sameall = OUTDAT.All_issame(indsthis);
                Targall = OUTDAT.All_istarg(indsthis);
                Yall = [];
                tbin = OUTDAT.All_FRsmooth_t(find(indsthis,1, 'first'));
                for j=indsthis'
                    
                    FRmeanAll = OUTDAT.AllMinusBase_FRmeanAll{j};
                    FRsemAll = OUTDAT.AllMinusBase_FRsemAll{j};
                    tbin = OUTDAT.AllMinusBase_tbinAll{j};
                    
                    Yall = [Yall FRmeanAll{epochtoplot}];
                end
                
                % ############################# PLOT 1 - individual syls
                ymean_all = mean(Yall,2);
                Yall_diff = Yall - ymean_all;
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                ylabel({[bname], [ename] ['sw' num2str(ss) ', neur' num2str(nn)]});
                title('deviation from glob mean');
                % -- diff
                pcol = 'r';
                plot(tbin, Yall_diff(:, Sameall==0 & Targall==0), 'Color', pcol);
                
                % -- same
                try
                    pcol = 'b';
                    plot(tbin, Yall_diff(:, Sameall==1 & Targall==0), 'Color', pcol, 'LineWidth', 2);
                catch err
                end
                
                % -- TARG
                pcol = 'k';
                plot(tbin, Yall_diff(:, Sameall==1 & Targall==1), 'Color', pcol, 'LineWidth', 2);
                
                lt_plot_zeroline;
                
                
                % ==== for each one, get absolute value deviation from
                % global mean
                ymean_all = mean(Yall,2);
                Yall = abs(Yall - ymean_all);
                
                
                % ############################# PLOT 2 - individual syls
                % [abs balue]
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                ylabel({[bname], [ename] ['sw' num2str(ss) ', neur' num2str(nn)]});
                title('abs dev from glob mean');
                
                % -- diff
                pcol = 'r';
                plot(tbin, Yall(:, Sameall==0 & Targall==0), 'Color', pcol);
                
                % -- same
                try
                    pcol = 'b';
                    plot(tbin, Yall(:, Sameall==1 & Targall==0), 'Color', pcol, 'LineWidth', 2);
                catch err
                end
                
                % -- TARG
                pcol = 'k';
                plot(tbin, Yall(:, Sameall==1 & Targall==1), 'Color', pcol, 'LineWidth', 2);
                
                lt_plot_zeroline;
                
                
                % ############################# PLOT 3 - mean across syls
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                ylabel({[bname], [ename] ['sw' num2str(ss) ', neur' num2str(nn)]});
                title('abs dev from glob mean, avg within syl type');
                
                % -- diff
                pcol = 'r';
                y = Yall(:, Sameall==0 & Targall==0);
                shadedErrorBar(tbin, mean(y,2),lt_sem(y'), {'Color', pcol}, 1);
                
                % -- same
                pcol = 'b';
                y = Yall(:, Sameall==1 & Targall==0);
                if size(y,2)>1
                    shadedErrorBar(tbin, mean(y,2),lt_sem(y'), {'Color', pcol}, 1);
                else
                    plot(tbin, mean(y,2), 'Color', pcol);
                end
                
                % -- TARG
                pcol = 'k';
                y = Yall(:, Sameall==1 & Targall==1);
                if size(y,2)>1
                    shadedErrorBar(tbin, mean(y,2),lt_sem(y'), {'Color', pcol}, 1);
                else
                    plot(tbin, mean(y,2), 'Color', pcol);
                end
                
                lt_plot_zeroline;
                
            end
        end
    end
end
end

