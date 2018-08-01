%% lt 5/24/18 - does direction of AFP bias change during learning?
% ACROSS EXPERIMENTS FOR THE SAME BIRD/SYLLABLE?
% NOTE: sumary data for PBS and MUSC (summary across expts) uses mean of
% day means, with perfectly matched days across PBS and MUSC.

%% LT 6/20/16
function lt_seq_dep_pitch_ACROSSBIRDS_Hamish2(SeqDepPitch_AcrossBirds, Params, ...
    plotsametype, useHandLab, plotPC, collectAllSyls, PCtimeWindowUsingWN)

%% HAND ENTERED time windows for pitch contours
% NOTE: if there is WN during training, then this would be based on
% abseline, so there is a chance will affect measure during WN training.
% PCtimewindows = {'pu53wh88', 'SeqDepPitchLMAN', [19 21], ...
%
% PCtimewindows = {'pu11wh87', 'SeqDepPitchLMAN', [0.019 0.043], ...
%     'pu11wh87', 'SeqDepPitchLMAN2', [0.018 0.045], ...
%     'pu11wh87', 'SeqDepPitchLMAN3', [0.019 0.037], ...
%     'pu11wh87', 'SeqDepPitchLMAN6', [0.019 0.044], ...
%     'gr41gr90', 'SeqDepPitchLMAN', [0.02 0.04], ...
%     'gr41gr90', 'SeqDepPitchLMAN2', [0.021 0.046], ...
%     'rd23gr89', 'SeqDepPitchLMAN', [0.018 0.045], ...
%     'rd23gr89',  'SeqDepPitchLMAN2', [0.019 0.029], ...
%     'rd28pu64',  'SeqDepPitchLMAN', [0.023 0.047], ...
%     'rd28pu64',  'SeqDepPitchLMAN2', [0.024 0.049], ...
%     'bk34bk68',  'SeqDepPitchLMAN', [0.018 0.041], ...
%     'bk34bk68',  'SeqDepPitchLMAN3', [0.017 0.038], ...
%     'wh4wh77',  'SeqDepPitchLMAN', [0.029 0.043], ...
%     'wh25pk77',  'SeqDepPitchLMAN', [0.027 0.037]}; % bird, expt, window (in sec),
%

PCtimewindows = {...
    'pu11wh87', 'SeqDepPitchLMAN', 'b', [0.021 0.043], ...
    'pu11wh87', 'SeqDepPitchLMAN', 'aB', [0.023 0.047], ...
    'pu11wh87', 'SeqDepPitchLMAN', 'a', [0.058 0.081], ...
    'pu11wh87', 'SeqDepPitchLMAN', 'c', [0.023 0.07], ...
    'pu11wh87', 'SeqDepPitchLMAN', 'd', [0.018 0.033], ...
    'pu11wh87', 'SeqDepPitchLMAN2', 'b', [0.019 0.038], ...
    'pu11wh87', 'SeqDepPitchLMAN2', 'aB', [0.019 0.041], ...
    'pu11wh87', 'SeqDepPitchLMAN2', 'a', [0.057 0.07], ...
    'pu11wh87', 'SeqDepPitchLMAN2', 'd', [0.017 0.031], ...
    'pu11wh87', 'SeqDepPitchLMAN2', 'c', [0.017 0.055], ...
    'pu11wh87', 'SeqDepPitchLMAN3', 'b', [0.021 0.037], ...
    'pu11wh87', 'SeqDepPitchLMAN3', 'aB', [0.018 0.04], ...
    'pu11wh87', 'SeqDepPitchLMAN3', 'a', [], ...
    'pu11wh87', 'SeqDepPitchLMAN3', 'c', [0.027 0.065], ...
    'pu11wh87', 'SeqDepPitchLMAN3', 'd', [0.017 0.027], ...
    'pu11wh87', 'SeqDepPitchLMAN6', 'b', [0.02 0.038], ...
    'pu11wh87', 'SeqDepPitchLMAN6', 'a', [0.055 0.077], ...
    'pu11wh87', 'SeqDepPitchLMAN6', 'c', [0.023 0.065], ...
    'pu11wh87', 'SeqDepPitchLMAN6', 'd', [0.018 0.033], ...
    'gr41gr90', 'SeqDepPitchLMAN', 'b', [0.021 0.04], ...
    'gr41gr90', 'SeqDepPitchLMAN', 'c', [0.025 0.075], ...
    'gr41gr90', 'SeqDepPitchLMAN', 'a', [0.063 0.082], ...
    'gr41gr90', 'SeqDepPitchLMAN2', 'b', [0.021 0.042], ...
    'gr41gr90', 'SeqDepPitchLMAN2', 'jbbacB', [0.022 0.039], ...
    'gr41gr90', 'SeqDepPitchLMAN2', 'c', [0.028 0.08], ...
    'gr41gr90', 'SeqDepPitchLMAN2', 'a', [0.062 0.082], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'dbB', [0.021 0.038], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'dB', [0.021 0.042], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'cB', [0.021 0.04], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'c', [0.018 0.066], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'd', [0.041 0.058], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'a', [0.047 0.065], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'h', [], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'b', [0.023 0.037], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'd', [0.039 0.045], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'a', [0.044 0.055], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'g', [0.036 0.052], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'c', [0.017 0.055], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'k', [0.019 0.029], ...
    'rd28pu64',  'SeqDepPitchLMAN', 'b', [0.025 0.044], ...
    'rd28pu64',  'SeqDepPitchLMAN', 'k', [], ...
    'rd28pu64',  'SeqDepPitchLMAN', 'd', [0.032 0.037], ...
    'rd28pu64',  'SeqDepPitchLMAN', 'a', [0.053 0.063], ...
    'rd28pu64',  'SeqDepPitchLMAN2', 'b', [0.024 0.045], ...
    'rd28pu64',  'SeqDepPitchLMAN2', 'd', [0.031 0.037], ...
    'rd28pu64',  'SeqDepPitchLMAN2', 'a', [0.051 0.062], ...
    'rd28pu64',  'SeqDepPitchLMAN2', 'k', [], ...
    'bk34bk68',  'SeqDepPitchLMAN', 'jjB', [0.026 0.046], ...
    'bk34bk68',  'SeqDepPitchLMAN', 'ljB', [0.026 0.046], ...
    'bk34bk68',  'SeqDepPitchLMAN', 'jjbB', [0.02 0.038], ...
    'bk34bk68',  'SeqDepPitchLMAN', 'ljbB', [0.02 0.038], ...
    'bk34bk68',  'SeqDepPitchLMAN', 'a', [], ... % empty means ignore (noisy)
    'bk34bk68',  'SeqDepPitchLMAN3', 'jjB', [0.024 0.043], ...
    'bk34bk68',  'SeqDepPitchLMAN3', 'ljB', [0.024 0.043], ...
    'bk34bk68',  'SeqDepPitchLMAN3', 'jjbB', [0.02 0.038], ...
    'bk34bk68',  'SeqDepPitchLMAN3', 'ljbB', [0.02 0.038], ...
    'bk34bk68',  'SeqDepPitchLMAN3', 'a', [], ... % empty means ignore (noisy)
    'wh4wh77',  'SeqDepPitchLMAN', 'b', [0.029 0.041], ...
    'wh4wh77',  'SeqDepPitchLMAN', 'n', [], ...
    'wh4wh77',  'SeqDepPitchLMAN', 'c', [0.038 0.07], ...
    'wh4wh77',  'SeqDepPitchLMAN', 'k', [0.017 0.023], ...
    'wh4wh77',  'SeqDepPitchLMAN', 'a', [0.045 0.057], ...
    'wh4wh77',  'SeqDepPitchLMAN', 'd', [], ...
    'wh25pk77',  'SeqDepPitchLMAN', 'b', [0.027 0.037], ...
    'wh25pk77',  'SeqDepPitchLMAN', 'hB', [0.017 0.03], ...
    'wh25pk77',  'SeqDepPitchLMAN', 'kB', [0.017 0.03], ...
    'wh25pk77',  'SeqDepPitchLMAN', 'h', [0.03 0.041], ...
    'wh25pk77',  'SeqDepPitchLMAN', 'a', [0.042 0.069], ...
    }; % bird, expt, window (in sec),
%     'rd23gr89', 'SeqDepPitchLMAN', 'h', [0.033 0.055], ... ALTERNATIVE
%     'wh4wh77',  'SeqDepPitchLMAN', 'n', [0.053 0.059], ... % ALTERNATIVE
%     'wh4wh77',  'SeqDepPitchLMAN', 'd', [0.027 0.035], ... ALTE$RNATIVE.


PCtimewindows_WN = {...
    'pu11wh87', 'SeqDepPitchLMAN', 'bccB', [0.022 0.027], ...
    'pu11wh87', 'SeqDepPitchLMAN2', 'abB', [0.016 0.022], ...
    'pu11wh87', 'SeqDepPitchLMAN3', 'dccB', [0.021 0.037], ...
    'pu11wh87', 'SeqDepPitchLMAN6', 'aB', [0.021 0.038], ...
    'gr41gr90', 'SeqDepPitchLMAN', 'jBba', [0.024 0.04], ...
    'gr41gr90', 'SeqDepPitchLMAN2', 'jbBa', [0.024 0.042], ...
    'rd23gr89', 'SeqDepPitchLMAN', 'dB', [0.024 0.042], ...
    'rd23gr89',  'SeqDepPitchLMAN2', 'cB', [0.023 0.037], ...
    'rd28pu64',  'SeqDepPitchLMAN', 'kjB', [0.026 0.044], ...
    'rd28pu64',  'SeqDepPitchLMAN2', 'jjB', [0.025 0.044], ...
    'bk34bk68',  'SeqDepPitchLMAN', 'ljbB', [0.021 0.038], ...
    'bk34bk68',  'SeqDepPitchLMAN3', 'jjbB', [0.021 0.036], ...
    'wh4wh77',  'SeqDepPitchLMAN', 'cbB', [0.029 0.041], ...
    'wh25pk77',  'SeqDepPitchLMAN', 'hbB', [0.027 0.037], ...
    }; % bird, expt, window (in sec),

% if collectAllSyls==1
%     % ------ becuase currently only contains windows for targ and s-type
%     disp('USING PREVIOUS PITCH CONTOUR WINDOWS...');
%     PCtimewindows = {};
% end

%% ================= which hand coded time windows to use?

if PCtimeWindowUsingWN==1
   % then only gets target syllables, 
    % using timw windows that were defined looking at boht base and WN
    % trials.
    PCtimewindows = PCtimewindows_WN;
    disp('NOTE: only getting target syls for pitch contour. if ok then continue...');
    pause 
end
    

%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

%%
useLearningRelLastBlineDay=0; % otherwise will use rel entire baseline
useWNday2=0; % otherwise will use 1. using 2 allows to get some expt without songs on day 1
takeDay2forExptWithNoDay1=1; % otherwise will throw out those experiments.


%% ========================= PLOT EACH EXPERIMENTS

All_learndir = [];
All_basebiaspval = [];
All_FF_BASE_PBS = [];
All_FF_WN_PBS = [];
All_FF_BASE_MUSC = [];
All_FF_WN_MUSC = [];
All_Birdnum = [];

All_WNdayrange = [];

All_CV_BASE_PBS = [];
All_CV_WN_PBS = [];
All_CV_BASE_MUSC = [];
All_CV_WN_MUSC = [];

exptcount = 0;
All_exptcounter = [];
All_Istarg = [];
All_Issame = [];


% ------------ PC STUFF
All_PitchCont_BASE_PBS = [];
All_PitchCont_BASE_MUSC = [];
All_PitchCont_WN_PBS = [];
All_PitchCont_WN_MUSC = [];

count = 0;
for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    %     figcount=1;
    %     subplotrows=4;
    %     subplotcols=2;
    %     fignums_alreadyused=[];
    %     hfigs=[];
    
    for ii=1:numexpts
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        % =============================== GET SET OF SYLS TO COLLECT
        SylsToCollect = {};
        if useHandLab==0
            sametypesyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        elseif useHandLab==1
            % ------ go thru all syls and keep ones that are same type, but
            % not target.
            allsyls = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            sametypesyls = {};
            for j=1:length(allsyls)
                syltmp = allsyls{j};
                sametmp = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syltmp).similar_to_targ_HandLab;
                targtmp = strcmp(syltmp, targsyl);
                if sametmp==1 & targtmp==0
                    sametypesyls = [sametypesyls syltmp];
                end
            end
            disp([targsyl ' ====== ' sametypesyls]);
        end
        if collectAllSyls==1
            % --- then collect targ, same, diff
            SylsToCollect = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        else
            % ---- then collect only targ, same
            SylsToCollect = [{targsyl} sametypesyls];
        end
        
        % ================================= WHICH syls to plot online
        if plotsametype==0
            SylsToPlot={targsyl};
        elseif plotsametype==1
            SylsToPlot = [{targsyl} sametypesyls];
        end
        
        
        %% ====== PLOT RAW DAT FOR THIS DAY TO COMPARE TO EXTRACTED STATS
        plotLarge=1;
        BirdToPlot=birdname;
        ExptToPlot=exptname;
        overlayMeans=1;
        UseSylColors=0; % 0 is black;
        flipsign=1; % plot pos, if neg
        use_std=0; % applies to mean plots only!! (std, instead of sem)
        plotRawFF=1; % 1=raw FF, 0= baseline subtracted (MUSC-MUSC)
        OverlayLMANStats=1; % plots mean shift in single target learning window (defined below)
        OverlayMUSC_days=[];
        plotRunningCV = 0;
        lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds,...
            Params, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, ...
            plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, ...
            OverlayMUSC_days, plotLarge, plotRunningCV)
        
        count = count+1;
        
        
        %% ############################## COLLECT STATS
        DatThis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning;
        
        if isfield(DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'final_extracted_window')
            
            %             SylsToCollect = [{targsyl}];
            exptcount = exptcount+1;
            
            for ss = 1:length(SylsToCollect)
                sylthis = SylsToCollect{ss};
                istarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).is_target;
                if useHandLab==0
                    issame = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).similar_to_targ;
                elseif useHandLab==1
                    issame = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).similar_to_targ_HandLab;
                end
                
                %% ============================================= MUSC
                % ============================= BASE
                tvals = DatThis.AllDays_PlotLearning.EpochData_MUSC.Baseline.(sylthis).Tvals_WithinTimeWindow;
                % ----------- 1) MEAN PITCH
                basedays = unique(floor(tvals));
                firstday = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.SeqFilter.FirstDay;
                basedayinds = lt_convert_EventTimes_to_RelTimes(firstday, basedays);
                basedayinds = basedayinds.JustDays_rel;
                
                doMUSC = 1;
                [ffmean_BASE_MUSC, ffCV_BASE_MUSC, ffvals_base_MUSC] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                    i, ii, basedayinds, sylthis, doMUSC);
                % ------ 2) PITCH CONTOUR
                [OUTSTRUCT] = ...
                    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub2(SeqDepPitch_AcrossBirds, ...
                    i, ii, sylthis, 'MUSC', basedayinds, PCtimewindows, PCtimeWindowUsingWN);
                All_PitchCont_BASE_MUSC = [All_PitchCont_BASE_MUSC; OUTSTRUCT];
                
                
                % ============================ DUR WN
                WNdays = DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds;
                % ------------ 1) MEAN PITCH
                doMUSC = 1;
                [ffmean_WN_MUSC, ffCV_WN_MUSC] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                    i, ii, WNdays, sylthis, doMUSC);
                % ------------ 2) PITCH CONTOUR
                [OUTSTRUCT] = ...
                    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub2(SeqDepPitch_AcrossBirds, ...
                    i, ii, sylthis, 'MUSC', WNdays, PCtimewindows, PCtimeWindowUsingWN);
                All_PitchCont_WN_MUSC = [All_PitchCont_WN_MUSC; OUTSTRUCT];
                
                
                
                %% ============================================= PBS
                % ============================ BASE
                % ----- 1) MEAN PITCH
                doMUSC = 0;
                [ffmean_BASE_PBS, ffCV_BASE_PBS, ffvals_base_PBS] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                    i, ii, basedayinds, sylthis, doMUSC);
                % ------ 2) PITCH CONTOUR
                [OUTSTRUCT] = ...
                    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub2(SeqDepPitch_AcrossBirds, ...
                    i, ii, sylthis, 'PBS', basedayinds, PCtimewindows, PCtimeWindowUsingWN);
                All_PitchCont_BASE_PBS = [All_PitchCont_BASE_PBS; OUTSTRUCT];
                
                
                
                % =========== DURING WN
                % ----------- 1) MEAN PITCH
                doMUSC = 0;
                [ffmean_WN_PBS, ffCV_WN_PBS] = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub(SeqDepPitch_AcrossBirds, ...
                    i, ii, WNdays, sylthis, doMUSC);
                % ------------ 2) PITCH CONTOUR
                [OUTSTRUCT] = ...
                    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub2(SeqDepPitch_AcrossBirds, ...
                    i, ii, sylthis, 'PBS', WNdays, PCtimewindows, PCtimeWindowUsingWN);
                All_PitchCont_WN_PBS = [All_PitchCont_WN_PBS; OUTSTRUCT];
                
                
                
                
                %% ====================== PLOT EXTRACTED VALUES ON LEARNING
                % --- PBS
                line([basedayinds(1)-0.5 basedayinds(end)+0.5], [ffmean_BASE_PBS ffmean_BASE_PBS], 'Color', 'b', 'LineWidth', 3);
                line([WNdays(1)-0.5 WNdays(end)+0.5], [ffmean_WN_PBS ffmean_WN_PBS], 'Color', 'b', 'LineWidth', 3);
                
                % --- MUSC
                line([basedayinds(1)-0.5 basedayinds(end)+0.5], [ffmean_BASE_MUSC ffmean_BASE_MUSC], 'Color', 'm', 'LineWidth', 3);
                line([WNdays(1)-0.5 WNdays(end)+0.5], [ffmean_WN_MUSC ffmean_WN_MUSC], 'Color', 'm', 'LineWidth', 3);
                
                
                %% ======================== PLOT PC?
                if plotPC ==1
                    Ntoplot = 20; % num renditions to contour;
                    if istarg==1
                        hsplots = [];
                        lt_figure; hold on;
                        ncol = 3;
                        
                        % ############################ BASELINE, PBS
                        pcol = 'k';
                        pcmat = All_PitchCont_BASE_PBS(end).All_PCmat{end};
                        twind = All_PitchCont_BASE_PBS(end).All_twind(end,:);
                        tbins = All_PitchCont_BASE_PBS(end).All_tbins{end};
                        ffmat = All_PitchCont_BASE_PBS(end).All_ffvals{end};
                        tt = tbins(twind);
                        
                        % ---------- 1) LOWER 20
                        hsplot = lt_subplot(4,ncol,1); hold on;
                        hsplots = [hsplots hsplot];
                        title('BASE, PBS, lower ff');
                        
                        [~, indsort] = sort(ffmat);
                        pcthis = pcmat(indsort(1:Ntoplot), :);
                        plot(tbins, pcthis', '-', 'Color', pcol);
                        % -- lines for twind
                        line([tt(1) tt(1)], ylim, 'Color', 'r');
                        line([tt(end) tt(end)], ylim, 'Color', 'r');
                        
                        % ---------- 1) HIGHER 20
                        hsplot = lt_subplot(4,ncol,2); hold on;
                        hsplots = [hsplots hsplot];
                        title('BASE, PBS, higher ff');
                        
                        [~, indsort] = sort(ffmat);
                        pcthis = pcmat(indsort(end-Ntoplot+1:end), :);
                        plot(tbins, pcthis', '-k', 'Color', pcol);
                        % -- lines for twind
                        line([tt(1) tt(1)], ylim, 'Color', 'r');
                        line([tt(end) tt(end)], ylim, 'Color', 'r');
                        
                        % ----------- OVERLAY MEANS
                        hsplot = lt_subplot(4, ncol, 3); hold on;
                        hsplots = [hsplots, hsplot];
                        title('overlaid means');
                        % - lower 20 (PBS)
                        pcthis = pcmat(indsort(1:Ntoplot), :);
                        pcmean = mean(pcthis,1);
                        pcsem = lt_sem(pcthis);
                        shadedErrorBar(tbins, pcmean, pcsem, {'Color', [0.6 0.6 0.6]}, 1);
                        
                        % - higher 20 (PBS)
                        pcthis = pcmat(indsort(end-Ntoplot+1:end), :);
                        pcmean = mean(pcthis,1);
                        pcsem = lt_sem(pcthis);
                        shadedErrorBar(tbins, pcmean, pcsem, {'Color', [0.6 0.6 0.6]}, 1);
                        
                        
                        % ############################ BASELINE, MUSC
                        pcol = 'r';
                        pcmat = All_PitchCont_BASE_MUSC(end).All_PCmat{end};
                        twind = All_PitchCont_BASE_MUSC(end).All_twind(end,:);
                        tbins = All_PitchCont_BASE_MUSC(end).All_tbins{end};
                        ffmat = All_PitchCont_BASE_MUSC(end).All_ffvals{end};
                        tt = tbins(twind);
                        
                        if size(pcmat,1)>Ntoplot
                            % ---------- 1) LOWER 20
                            hsplot = lt_subplot(4,ncol,4); hold on;
                            hsplots = [hsplots hsplot];
                            title('BASE, lower ff');
                            
                            [~, indsort] = sort(ffmat);
                            pcthis = pcmat(indsort(1:Ntoplot), :);
                            plot(tbins, pcthis', '-', 'Color', pcol);
                            % -- lines for twind
                            line([tt(1) tt(1)], ylim, 'Color', 'r');
                            line([tt(end) tt(end)], ylim, 'Color', 'r');
                            
                            % ---------- 1) HIGHER 20
                            hsplot = lt_subplot(4,ncol,5); hold on;
                            hsplots = [hsplots hsplot];
                            title('BASE, higher ff');
                            
                            [~, indsort] = sort(ffmat);
                            pcthis = pcmat(indsort(end-Ntoplot+1:end), :);
                            plot(tbins, pcthis', '-k', 'Color', pcol);
                            % -- lines for twind
                            line([tt(1) tt(1)], ylim, 'Color', 'r');
                            line([tt(end) tt(end)], ylim, 'Color', 'r');
                        end
                        
                        % - OVERLAY MEANS - MUSC (all)
                        lt_subplot(4, ncol, 3); hold on;
                        pcthis = pcmat;
                        pcmean = mean(pcthis,1);
                        pcsem = lt_sem(pcthis);
                        shadedErrorBar(tbins, pcmean, pcsem, {'Color', 'r'}, 1);
                        
                        
                        
                        % ############################ TRAIN, PBS
                        pcol = 'k';
                        pcmat = All_PitchCont_WN_PBS(end).All_PCmat{end};
                        twind = All_PitchCont_WN_PBS(end).All_twind(end,:);
                        tbins = All_PitchCont_WN_PBS(end).All_tbins{end};
                        ffmat = All_PitchCont_WN_PBS(end).All_ffvals{end};
                        tt = tbins(twind);
                        
                        if size(pcmat,1)>Ntoplot
                            
                            % ---------- 1) LOWER 20
                            hsplot = lt_subplot(4,ncol,7); hold on;
                            hsplots = [hsplots hsplot];
                            title('WN, lower ff');
                            
                            [~, indsort] = sort(ffmat);
                            pcthis = pcmat(indsort(1:Ntoplot), :);
                            plot(tbins, pcthis', '-', 'Color', pcol);
                            % -- lines for twind
                            line([tt(1) tt(1)], ylim, 'Color', 'r');
                            line([tt(end) tt(end)], ylim, 'Color', 'r');
                            
                            % ---------- 1) HIGHER 20
                            hsplot = lt_subplot(4,ncol,8); hold on;
                            hsplots = [hsplots hsplot];
                            title('WN, higher ff');
                            
                            [~, indsort] = sort(ffmat);
                            pcthis = pcmat(indsort(end-Ntoplot+1:end), :);
                            plot(tbins, pcthis', '-k', 'Color', pcol);
                            % -- lines for twind
                            line([tt(1) tt(1)], ylim, 'Color', 'r');
                            line([tt(end) tt(end)], ylim, 'Color', 'r');
                            
                            % ----------- OVERLAY MEANS
                            hsplot = lt_subplot(4, ncol, 9); hold on;
                            hsplots = [hsplots, hsplot];
                            title('overlaid means');
                            % - lower 20 (PBS)
                            pcthis = pcmat(indsort(1:Ntoplot), :);
                            pcmean = mean(pcthis,1);
                            pcsem = lt_sem(pcthis);
                            shadedErrorBar(tbins, pcmean, pcsem, {'Color', [0.6 0.6 0.6]}, 1);
                            
                            % - higher 20 (PBS)
                            pcthis = pcmat(indsort(end-Ntoplot+1:end), :);
                            pcmean = mean(pcthis,1);
                            pcsem = lt_sem(pcthis);
                            shadedErrorBar(tbins, pcmean, pcsem, {'Color', [0.6 0.6 0.6]}, 1);
                            
                            
                        end
                        
                        % ############################ TRAIN, MUSC
                        pcol = 'r';
                        pcmat = All_PitchCont_WN_MUSC(end).All_PCmat{end};
                        twind = All_PitchCont_WN_MUSC(end).All_twind(end,:);
                        tbins = All_PitchCont_WN_MUSC(end).All_tbins{end};
                        ffmat = All_PitchCont_WN_MUSC(end).All_ffvals{end};
                        tt = tbins(twind);
                        
                        if size(pcmat,1)>Ntoplot
                            
                            % ---------- 1) LOWER 20
                            hsplot = lt_subplot(4,ncol,10); hold on;
                            hsplots = [hsplots hsplot];
                            title('WN, lower ff');
                            
                            [~, indsort] = sort(ffmat);
                            pcthis = pcmat(indsort(1:Ntoplot), :);
                            plot(tbins, pcthis', '-', 'Color', pcol);
                            % -- lines for twind
                            line([tt(1) tt(1)], ylim, 'Color', 'r');
                            line([tt(end) tt(end)], ylim, 'Color', 'r');
                            
                            % ---------- 1) HIGHER 20
                            hsplot = lt_subplot(4,ncol,11); hold on;
                            hsplots = [hsplots hsplot];
                            title('WN, higher ff');
                            
                            [~, indsort] = sort(ffmat);
                            pcthis = pcmat(indsort(end-Ntoplot+1:end), :);
                            plot(tbins, pcthis', '-k', 'Color', pcol);
                            % -- lines for twind
                            line([tt(1) tt(1)], ylim, 'Color', 'r');
                            line([tt(end) tt(end)], ylim, 'Color', 'r');
                            
                        end
                        
                        % - OVERLAY MEANS - MUSC (all)
                        lt_subplot(4, ncol, 9); hold on;
                        pcthis = pcmat;
                        pcmean = mean(pcthis,1);
                        pcsem = lt_sem(pcthis);
                        shadedErrorBar(tbins, pcmean, pcsem, {'Color', 'r'}, 1);
                        
                        
                        
                        % ====================== format
                        linkaxes(hsplots, 'xy');
                        
                        
                    end
                end
                
                %% ############################# COLLECT STATS
                
                learndir = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                All_learndir = [All_learndir; learndir];
                
                [~, p] = ttest2(ffvals_base_MUSC, ffvals_base_PBS);
                if p<0.05
                    lt_plot_pvalue(p, 'base bias (ttest)', 1);
                end
                All_basebiaspval = [All_basebiaspval; p];
                
                All_FF_BASE_PBS = [All_FF_BASE_PBS; ffmean_BASE_PBS];
                All_FF_WN_PBS = [All_FF_WN_PBS; ffmean_WN_PBS];
                All_FF_BASE_MUSC = [All_FF_BASE_MUSC; ffmean_BASE_MUSC];
                All_FF_WN_MUSC = [All_FF_WN_MUSC; ffmean_WN_MUSC];
                
                All_CV_BASE_PBS = [All_CV_BASE_PBS; ffCV_BASE_PBS];
                All_CV_WN_PBS = [All_CV_WN_PBS; ffCV_WN_PBS];
                All_CV_BASE_MUSC = [All_CV_BASE_MUSC; ffCV_BASE_MUSC];
                All_CV_WN_MUSC = [All_CV_WN_MUSC; ffCV_WN_MUSC];
                
                
                All_Birdnum = [All_Birdnum; i];
                All_exptcounter = [All_exptcounter; exptcount];
                
                WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd + ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
                All_WNdayrange = [All_WNdayrange; [WNdays(1) WNdays(end)]-WNday1+1];
                
                All_Istarg = [All_Istarg; istarg];
                All_Issame = [All_Issame; issame];
            end
            
            
        else
            lt_plot_annotation(1, 'NOT COLLECTING DATA! no dat in final window...', 'm');
        end
        
        
    end
end

disp(['Num expts plotted: ' num2str(count)]);


%% ################### PUT ALL DATA INTO STRUCTURE

DATSTRUCT.All_Birdnum = All_Birdnum;
DATSTRUCT.All_CV_BASE_MUSC = All_CV_BASE_MUSC;
DATSTRUCT.All_CV_BASE_PBS = All_CV_BASE_PBS;
DATSTRUCT.All_CV_WN_MUSC = All_CV_WN_MUSC;
DATSTRUCT.All_CV_WN_PBS = All_CV_WN_PBS;
DATSTRUCT.All_FF_BASE_MUSC = All_FF_BASE_MUSC;
DATSTRUCT.All_FF_BASE_PBS = All_FF_BASE_PBS;
DATSTRUCT.All_FF_WN_MUSC = All_FF_WN_MUSC;
DATSTRUCT.All_FF_WN_PBS = All_FF_WN_PBS;
DATSTRUCT.All_Istarg = All_Istarg;
DATSTRUCT.All_Issame = All_Issame;
DATSTRUCT.All_PitchCont_BASE_MUSC = All_PitchCont_BASE_MUSC;
DATSTRUCT.All_PitchCont_BASE_PBS = All_PitchCont_BASE_PBS;
DATSTRUCT.All_PitchCont_WN_MUSC = All_PitchCont_WN_MUSC;
DATSTRUCT.All_PitchCont_WN_PBS = All_PitchCont_WN_PBS;
DATSTRUCT.All_WNdayrange = All_WNdayrange;
DATSTRUCT.All_basebiaspval = All_basebiaspval;
DATSTRUCT.All_exptcounter = All_exptcounter;
DATSTRUCT.All_learndir = All_learndir;


%% ########## COLLECT WIGGLE AND BASELINE BIAS [recalculated using new time window]
% NOTE: if multiple days, then days mean across days of day stats.
% NOTE: calcualte baseline bias using pitch countours and update time
% windows

% Indstoplot = find(All_Istarg==1);
Indstoplot = find(All_Issame==1 | All_Issame==0); % uses all syllables
NrendsORIG = 25; % renditions to take from edges (will minimize if not enough trials
usemedian = 1; % if 0, then uses means. if 1, then uses MAD and medians


% ============== VERSION 1 - using wiggle normalzied to MUSC wiggle
normwiggle = 1; % if 1, then normalized wiggle to the wiggle for MUSC (matched sample size)
DATSTRUCT = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub7(DATSTRUCT, Indstoplot, ...
    NrendsORIG, usemedian, normwiggle);


% ============== VERSION 2 - 
normwiggle = 0; 
DATSTRUCT = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub7(DATSTRUCT, Indstoplot, ...
    NrendsORIG, usemedian, normwiggle);


% ################################### IF WANT TO ALSO INCLUDE SAME
% MESASUREMNTS BUT DURING WN TRAINING ...
normwiggle = 0; % set to 0, since want to compare WN and base, so don't want 
% to normalize to MUSC

% FIRST, EXTRACT BASELINE
DATSTRUCT = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub7(DATSTRUCT, Indstoplot, ...
    NrendsORIG, usemedian, normwiggle);

% SECOND, EXTRACT DURING TRAINING, AND ADD RELEVANT FIELDS TO DATSTRUCT
doBase = 0;
DATSTRUCT_TMP = lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub7(DATSTRUCT, Indstoplot, ...
    NrendsORIG, usemedian, normwiggle, doBase);

% THIRD, ADD RELEVANT FIELDS FROM TRAINING TO BASELINE
DATSTRUCT.Wiggle_WN = DATSTRUCT_TMP.Wiggle;


%% ######################################## PITCH CONTOUR STUFF [WIGGLE]
% ====== INCLUDES ALL SYLS (so can potentailly pseudoreplicate if have
% multiple expriments in same bird - will need to account for that)]

% NOTE: can try running twice based on the two versions above.
lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub6(DATSTRUCT);



%% ================= [PITCH CONTOUR/WIGGLE] SANITY CHECLS
% Go into script and run;
% NOTE: 
% section 1: useful (plots raw contours for all syls)
% section 2 (PBS vs MUSC) and section 3 (PBS vs PBS) less useful, since
% subseumed by RAW plots below...


if (0)
lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub8;
end

%% ========================= [PLOT RAW DATA - CHOOSE SPECIFIC SYLLABLE]
figcount=1;
subplotrows=4;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


Nrends = NrendsORIG;
indthis = 1;
dd =1; % which day

% ############################################ COMPARED TO MUSCIMOL
default_type = 'MUSC';
[fignums_alreadyused, hfigs, figcount, hsplot] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub5(default_type, DATSTRUCT, ...
    indthis, Nrends, dd, ...
    subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount)


% ############################################ COMPARED TO PBS
default_type = 'PBS';
[fignums_alreadyused, hfigs, figcount, hsplot] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub5(default_type, DATSTRUCT, ...
    indthis, Nrends, dd, ...
    subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount)


%% ========================= [PLOT RAW DATA FOR DIFFERENT SYLLABLES AS EXAMPLES]
% NOTE: same as above, except iterate over syllables of choice.

% ============================ MOIDIFY
plotrends = 'low';
% high or low, plots all syls, starting from eitehr high or low...

% ====================== PLOT ALL, WITH PAUSES IN BETWEEN
biasall = DATSTRUCT.BaseBiasAll;
[~, IndSort] = sort(biasall);

if strcmp(plotrends, 'low')
indstoplot = IndSort';
elseif strcmp(plotrends, 'high')
% --- go from high to lo
indstoplot = fliplr(indstoplot);
end

% =================================== RUNS
figcount=1;
subplotrows=4;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

dd =1; % which day
Nrends = NrendsORIG;

count = 0;
for i = indstoplot
    indthis = i;
    
    if any(isnan(DATSTRUCT.All_PitchCont_BASE_PBS(indthis).All_twind(1,:)))
        continue
    end
    
    % ############################################ COMPARED TO MUSCIMOL
    default_type = 'MUSC';
    [fignums_alreadyused, hfigs, figcount, hsplot] = ...
        lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub5(default_type, DATSTRUCT, ...
        indthis, Nrends, dd, ...
        subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    
    
    % ############################################ COMPARED TO PBS
    default_type = 'PBS';
    [fignums_alreadyused, hfigs, figcount, hsplot] = ...
        lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub5(default_type, DATSTRUCT, ...
        indthis, Nrends, dd, ...
        subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    
    % ------------ DO PAUSE?
    count = count+1;
    if mod(count, 8)==0
        pause
        close all;
    end
    
    
end

%% ==================================== is ff distribution skewed in direction of bias?

biasall =[];
skewall = [];
for j=1:length(DATSTRUCT.All_Birdnum)
    
    ff = DATSTRUCT.All_PitchCont_BASE_PBS(j).All_ffvals{1};
    ffmusc = DATSTRUCT.All_PitchCont_BASE_MUSC(j).All_ffvals{1};
    
    % --- get bias
    afpbias = mean(ff) - mean(ffmusc);
    
    
    % ---- get metric of skew
    if (1)
        p = prctile(ff, [5 50 95]);
    skew = (p(3) - p(2)) - (p(2) - p(1));
    else
        
    end
    
    % -------------------- COOLECT
    biasall = [biasall; afpbias];
    skewall = [skewall; skew];
    
end


lt_figure; hold on;
xlabel('baseline bias');
ylabel('skew');
lt_regress(skewall, biasall, 1);



%% ############################ TARGET SYL ONLY

lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub3(DATSTRUCT)




%% ########################### TARG + SAMETYPE

lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub4(DATSTRUCT)










