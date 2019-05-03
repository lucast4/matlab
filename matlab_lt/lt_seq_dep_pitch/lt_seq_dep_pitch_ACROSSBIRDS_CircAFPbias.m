%% lt 5/24/18 - does direction of AFP bias change during learning?
% ACROSS EXPERIMENTS FOR THE SAME BIRD/SYLLABLE?
% NOTE: sumary data for PBS and MUSC (summary across expts) uses mean of
% day means, with perfectly matched days across PBS and MUSC.

%% LT 6/20/16
function lt_seq_dep_pitch_ACROSSBIRDS_CircAFPbias(SeqDepPitch_AcrossBirds, Params, ...
    useHandLab)

%% ================= which hand coded time windows to use?

%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


%% ========================= PLOT EACH EXPERIMENTS

% ================
DATSTRUCT.birdnum = [];
DATSTRUCT.enum = [];
DATSTRUCT.syl = {};
DATSTRUCT.istarg = [];
DATSTRUCT.issame = [];
DATSTRUCT.learndirTarg = [];
DATSTRUCT.tvalsDays_PBSMUSC = {};
DATSTRUCT.ffvalsDays_PBSMUSC = {};
DATSTRUCT.ffMean_BaseWN_PbsMusc = {};

DATSTRUCT.BaseBiasPval = [];

DATSTRUCT.WNepochDays = {};
DATSTRUCT.BaseDays = {};
DATSTRUCT.sylnum = [];

for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname = SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
        ename = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
%         % =============================== GET SET OF SYLS TO COLLECT
%         if useHandLab==0
%             sametypesyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
%                 SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
%         elseif useHandLab==1
%             % ------ go thru all syls and keep ones that are same type, but
%             % not target.
%             allsyls = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%             
%             sametypesyls = {};
%             for j=1:length(allsyls)
%                 syltmp = allsyls{j};
%                 sametmp = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syltmp).similar_to_targ_HandLab;
%                 targtmp = strcmp(syltmp, targsyl);
%                 if sametmp==1 & targtmp==0
%                     sametypesyls = [sametypesyls syltmp];
%                 end
%             end
%             disp([targsyl ' ====== ' sametypesyls]);
%         end
        
        % --- then collect targ, same, diff
        SylsToCollect = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        figcount=1;
        subplotrows=3;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];


        %% ############################## COLLECT STATS
        DatThis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning;
        
        if ~isfield(DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'final_extracted_window')
            disp('SKIPPING - no final_extracted_window extracted - I assume it is bad expt');
            continue           
        end
        
        for ss = 1:length(SylsToCollect)
            
            
            
            % ================= METADAT
            sylthis = SylsToCollect{ss};
            istarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).is_target;
            if useHandLab==0
                issame = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).similar_to_targ;
            elseif useHandLab==1
                issame = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).similar_to_targ_HandLab;
            end
            
            
            learndirTarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            firstday=DatThis.Params.SeqFilter.FirstDay;
            
            
            %% ========================= RAW PITCH OVER ALL DAYS
            NumDays=length(DatThis.AllDays_PlotLearning.DataMatrix.(sylthis).Tvals);
            
            tvalsDays_PBSMUSC = cell(NumDays, 2); % days x [PBS, MUSC]
            ffvalsDays_PBSMUSC = cell(NumDays, 2); % days x [PBS, MUSC]
            
            for day=1:NumDays
                % ################################### PBS
                % then get all data from this day, not just time window
                tvals=cell2mat(DatThis.AllDays_PlotLearning.DataMatrix.(sylthis).Tvals{day});
                ffvals=cell2mat(DatThis.AllDays_PlotLearning.DataMatrix.(sylthis).FFvals{day});
                
                if isempty(ffvals)
                    assert(isempty(tvals))
                    continue
                end
                
                % --- convert tvals to day
                tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                tvals=tvals.FinalValue;
                
                % ======== OUTPUT
                tvalsDays_PBSMUSC{day, 1} = tvals; % days x [PBS, MUSC]
                ffvalsDays_PBSMUSC{day, 1} = ffvals; % days x [PBS, MUSC]
                
                % ################################### MUSC
                % then get all data from this day, not just time window
                tvals=cell2mat(DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis).Tvals{day});
                ffvals=cell2mat(DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(sylthis).FFvals{day});
                
                if isempty(ffvals)
                    assert(isempty(tvals))
                    continue
                end
                
                % - convert tvals to day
                tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                tvals=tvals.FinalValue;
                
                % ======== OUTPUT
                tvalsDays_PBSMUSC{day, 2} = tvals; % days x [PBS, MUSC]
                ffvalsDays_PBSMUSC{day, 2} = ffvals; % days x [PBS, MUSC]
                
                
                % ########################## MUSC SCHEDULE
                %                                 % ===== plot lines for MUSC switch + lag
                %                 lagtime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.Lag_time;
                %                 muscSchedule=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds{day};
                %
                %                 muscStart=day+muscSchedule.start/24;
                %                 muscEnd=day+muscSchedule.end/24;
                
            end
            
            
            %% ======================== AFP BIAS BASELINE AND WN
            
            ffMean_BaseWN_PbsMusc = nan(2,2);
            
            % ###################################################### MUSC
            firstday = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.SeqFilter.FirstDay;
            
            % ============================= BASE
            tvals = DatThis.AllDays_PlotLearning.EpochData_MUSC.Baseline.(sylthis).Tvals_WithinTimeWindow;
            % ----------- 1) MEAN PITCH
            basedays = unique(floor(tvals));
            basedayinds = lt_convert_EventTimes_to_RelTimes(firstday, basedays);
            basedayinds = basedayinds.JustDays_rel;
            
            doMUSC = 1;
            [ffmean_BASE_MUSC, ffCV_BASE_MUSC, ffvals_base_MUSC] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                i, ii, basedayinds, sylthis, doMUSC);
%             % ------ 2) PITCH CONTOUR
%             [OUTSTRUCT] = ...
%                 lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub2(SeqDepPitch_AcrossBirds, ...
%                 i, ii, sylthis, 'MUSC', basedayinds, PCtimewindows, PCtimeWindowUsingWN);
%             All_PitchCont_BASE_MUSC = [All_PitchCont_BASE_MUSC; OUTSTRUCT];
            ffMean_BaseWN_PbsMusc(1, 2) = ffmean_BASE_MUSC;
            
            
            % ============================ DUR WN
            WNdays = DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds;
            % ------------ 1) MEAN PITCH
            doMUSC = 1;
            [ffmean_WN_MUSC, ffCV_WN_MUSC] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                i, ii, WNdays, sylthis, doMUSC);
            ffMean_BaseWN_PbsMusc(2, 2) = ffmean_WN_MUSC;
            
            
            % ###################################################### PBS
            % ============================ BASE
            % ----- 1) MEAN PITCH
            doMUSC = 0;
            [ffmean_BASE_PBS, ffCV_BASE_PBS, ffvals_base_PBS] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                i, ii, basedayinds, sylthis, doMUSC);
            ffMean_BaseWN_PbsMusc(1, 1) = ffmean_BASE_PBS;
            

            % =========== DURING WN
            % ----------- 1) MEAN PITCH
            doMUSC = 0;
            [ffmean_WN_PBS, ffCV_WN_PBS] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                i, ii, WNdays, sylthis, doMUSC);
            ffMean_BaseWN_PbsMusc(2, 1) = ffmean_WN_PBS;
            
            
            
            
            %% ====================== PLOT EXTRACTED VALUES ON LEARNING
%             if istarg==1
%                 % --- PBS
%                 line([basedayinds(1)-0.5 basedayinds(end)+0.5], [ffmean_BASE_PBS ffmean_BASE_PBS], 'Color', 'b', 'LineWidth', 3);
%                 line([WNdays(1)-0.5 WNdays(end)+0.5], [ffmean_WN_PBS ffmean_WN_PBS], 'Color', 'b', 'LineWidth', 3);
%                 
%                 % --- MUSC
%                 line([basedayinds(1)-0.5 basedayinds(end)+0.5], [ffmean_BASE_MUSC ffmean_BASE_MUSC], 'Color', 'm', 'LineWidth', 3);
%                 line([WNdays(1)-0.5 WNdays(end)+0.5], [ffmean_WN_MUSC ffmean_WN_MUSC], 'Color', 'm', 'LineWidth', 3);
%             end

            %% ############################# COLLECT STATS
            
            [~, p] = ttest2(ffvals_base_MUSC, ffvals_base_PBS);
            DATSTRUCT.BaseBiasPval = [DATSTRUCT.BaseBiasPval; p];
            
            DATSTRUCT.birdnum = [DATSTRUCT.birdnum; i];
            DATSTRUCT.enum = [DATSTRUCT.enum; ii];
            DATSTRUCT.syl = [DATSTRUCT.syl; sylthis];
            DATSTRUCT.sylnum = [DATSTRUCT.sylnum; ss];
            DATSTRUCT.istarg = [DATSTRUCT.istarg; istarg];
            DATSTRUCT.issame = [DATSTRUCT.issame; issame];
            DATSTRUCT.learndirTarg = [DATSTRUCT.learndirTarg; learndirTarg];
            DATSTRUCT.tvalsDays_PBSMUSC = [DATSTRUCT.tvalsDays_PBSMUSC; tvalsDays_PBSMUSC];
            DATSTRUCT.ffvalsDays_PBSMUSC = [DATSTRUCT.ffvalsDays_PBSMUSC; ffvalsDays_PBSMUSC];
            DATSTRUCT.ffMean_BaseWN_PbsMusc = [DATSTRUCT.ffMean_BaseWN_PbsMusc; ffMean_BaseWN_PbsMusc];

            WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd + ...
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
%             DATSTRUCT.WNepochDays = [DATSTRUCT.WNepochDays; [WNdays(1) WNdays(end)]-WNday1+1];
            
            tvals = DatThis.AllDays_PlotLearning.EpochData.Baseline.(sylthis).Tvals_WithinTimeWindow;
            basedays = unique(floor(tvals));
            basedayinds = lt_convert_EventTimes_to_RelTimes(firstday, basedays);
            basedayinds = basedayinds.JustDays_rel;
            
            DATSTRUCT.BaseDays = [DATSTRUCT.BaseDays; basedayinds];
%             assert(WNday1 == basedayinds(end)+1);

            
            %% =============================== PLOT EACH SYLLABEL
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            ylabel([birdname '-' ename]);
            title([sylthis]);
            
            % -- PBS
            tvals = tvalsDays_PBSMUSC(basedayinds, 1);
            ffvals = ffvalsDays_PBSMUSC(basedayinds, 1);
            tvals = cell2mat(cellfun(@(x)x', tvals, 'UniformOutput', 0));
            ffvals = cell2mat(cellfun(@(x)x', ffvals, 'UniformOutput', 0));
            
            plot(tvals, ffvals, 'ok');
            
            % ------- PLOT SMOOTHED LINES FOR EACH DAY
            tvals = tvalsDays_PBSMUSC(basedayinds, 1);
            ffvals = ffvalsDays_PBSMUSC(basedayinds, 1);
            nday = length(tvals);
            for j=1:nday
               if length(ffvals{j})>(15+10)
                   t = lt_running_stats(tvals{j}, 15);
               f = lt_running_stats(ffvals{j}, 15);
                              shadedErrorBar(t.Mean, f.Mean, f.SEM, {'Color', 'k'}, 1);
               end
            end
            
            
            % -- MUSC
            tvals = tvalsDays_PBSMUSC(basedayinds, 2);
            ffvals = ffvalsDays_PBSMUSC(basedayinds, 2);
            tvals = cell2mat(cellfun(@(x)x', tvals, 'UniformOutput', 0));
            ffvals = cell2mat(cellfun(@(x)x', ffvals, 'UniformOutput', 0));
            
            plot(tvals, ffvals, 'or');
            
            
                        % --------------- OVERLAY DIRECTION OF BIAS
            y = [ffMean_BaseWN_PbsMusc(1, 1) ffMean_BaseWN_PbsMusc(1, 2)];
            x = [max(tvals)+0.1 max(tvals)+0.2];
            plot(x,y, '-')
            lt_plot(x(1), y(1), {'Color', 'k'});
            lt_plot(x(2), y(2), {'Color', 'r'});
            
            
            
        end
        
        
        
    end
end


%% ############ PLOT - each syllable




