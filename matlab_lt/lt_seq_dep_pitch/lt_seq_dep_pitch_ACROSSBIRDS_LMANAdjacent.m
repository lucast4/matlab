function lt_seq_dep_pitch_ACROSSBIRDS_LMANAdjacent(SeqDepPitch_AcrossBirds, PARAMS, ...
    norm_by_targsyl, epochfield_input, UseBaselineForCV, DispEachSylCVpval)
%% PARAMS
epochfield=epochfield_input;

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var')
    norm_by_targsyl=1;
end


%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

disp('--');
disp('Removing syllables that should not be analyzed - i.e. WN overlap, since not using catch songs. REMOVED:');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(PARAMS.global.LMAN.SylsToRemove_SingleDir, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+1};
            syls_toremove=PARAMS.global.LMAN.SylsToRemove_SingleDir{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
                    syls_unique(ind_to_remove)=[];
                    
                    % Put back into structure
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=syls_unique;
                    
                    % tell user this is done
                    disp([birdname '-' exptname ': removed ' tmp_sylremove]);
                end
            end
        end
    end
end


%% === list names of all types of syls
if (0)
    disp(' --- ');
    
    for i=1:length(SeqDepPitch_AcrossBirds.birds);
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        
        for ii=1:numexpts;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                    continue
                end
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1 ...
                        & SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl==0;
                    
                    
                    disp([birdname '-' exptname '-' syl]);
                end
            end
            
            
            
        end
    end
end

%% ==== list birds and expts
disp(' =============== birds/expts');
for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        disp([birdname '-' exptname]);
        
    end
end


%% [ COLLECT DATA] - PLOT AFP BIAS VS. MP LEARNING, DISTRIBUTIONS, AND OTHER THINGS
% lt_figure;
% hold on;
% title('
% line([0 0.05], [0 0.05]);

% epochfield=epochfield_input;

LearningPBS_all=[];

MPbias_all=[];
AFPbias_all=[];
SimDiff_all=[];
TargStatus_all=[];
PreSylSimilar_all=[];
Expt_count_all=[];
Yexpt_all={};
Ybird_all={};
Y_PosRelTarg_All=[];

Generalization_MP_all=[];
Generalization_AFP_all=[];
Generalization_Learn_all=[];

cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[];
cvRatio_pvalue_UsingAllVals_ALLEXPTS=[];

CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[];
CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[];

cvPBS_alldays_ALLEXPTS={};
cvMUSC_alldays_ALLEXPTS={};

TargLearnDirAll=[];


expt_count=1;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        
        % ==== ONE FIGURE PER EXPT
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % ========== COLLECT DATA
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
            
            disp(['COLLLECTED DATA FOR : ' birdname '-' exptname]);
            
            Y_FFmean_pbs=[];
            Y_FFmean_musc=[];
            Y_FFsem_pbs=[];
            Y_FFsem_musc=[];
            Y_syls={};
            Y_similar_diff=[];
            Y_istarg=[];
            Y_AFP_bias=[];
            Y_AcousticDist=[];
            Y_Corr=[];
            Y_presimilar=[];
            Yexpt={};
            Ybird={};
            Y_PosRelTarg=[];
            
            Y_Generalization_MP=[];
            Y_Generalization_AFP=[];
            Y_Generalization_Learn=[];
            
            
            % -- for CV stuff
            cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS=[];
            cvRatio_pvalue_UsingAllVals_ALLSYLS=[];
            
            CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS=[];
            CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS=[];
            
            cvPBS_alldays_ALLSYLS={};
            cvMUSC_alldays_ALLSYLS={};
            
            % ===== STATS AT TARGET
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            FF_PBS_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_pbs;
            FF_MUSC_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(targsyl).meanFF_musc;
            FF_AFP_targ=FF_PBS_targ-FF_MUSC_targ;
            
            targlearndir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                % MUSC)
                FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                
                FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                
                
                % ====== CALCULATE AFP AND MUSC RELATIVE TO
                % THOSE VALUES AT TARGET (I.E. GENERALIZATIONS)
                FF_AFP=FF_PBS-FF_MUSC;
                
                Generalization_MP=FF_MUSC/FF_MUSC_targ;
                Generalization_AFP=FF_AFP/FF_AFP_targ;
                Generalization_Learn=FF_PBS/FF_PBS_targ;
                
                
                
                % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                % ======================== calculate CV PBS and MUSC (for each day of
                if UseBaselineForCV==1;
                    TvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).Tvals_WithinTimeWindow;
                    TvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).Tvals_WithinTimeWindow;
                    
                    FFvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow;
                    FFvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).rawFF_WithinTimeWindow;
                    
                else
                    % inactivation, and mean across days)
                    TvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).TvalsWithinWindow_PBS;
                    TvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).TvalsWithinWindow_MUSC;
                    
                    FFvalsPBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).FFminusBaseWithinWindow_PBS_vals;
                    FFvalsMUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).FFminusBaseWithinWindow_MUSC_vals;
                    % convert FFvals to actual values, not diff from base
                    FFvalsPBS=FFvalsPBS+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                    FFvalsMUSC=FFvalsMUSC+SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                end
                
                % make sure tvals correspond to ffvals
                assert(length(TvalsPBS)==length(FFvalsPBS) & length(TvalsMUSC)==length(FFvalsMUSC), 'PRoblem - tvals and ffvals dont match');
                
                
                % -- for each day, calculate CV
                TvalsPBS_round=floor(TvalsPBS);
                TvalsMUSC_round=floor(TvalsMUSC);
                
                days=unique(TvalsPBS_round);
                
                cvPBS_alldays=[];
                cvMUSC_alldays=[];
                ffvals_DivideDayMean_AllDays_PBS=[];
                ffvals_DivideDayMean_AllDays_MUSC=[];
                % ffvals_MinusDayMean_AllDays_PBS=[];
                % ffvals_MinusDayMean_AllDays_MUSC=[];
                
                if isempty(days)
                    continue
                end
                
                for k=days
                    % for each day, get PBS and MUSC cv
                    
                    % --- PBS
                    inds=TvalsPBS_round==k; % only this day
                    ffvals=FFvalsPBS(inds);
                    cvPBS=std(ffvals)/abs(mean(ffvals));
                    
                    %                    ffvals_MinusMean=ffvals-mean(ffvals);
                    %                    ffvals_MinusDayMean_AllDays_PBS=[ffvals_MinusDayMean_AllDays_PBS ffvals_MinusMean];
                    
                    ffvals_DivideDayMean_AllDays_PBS=[ffvals_DivideDayMean_AllDays_PBS ffvals/mean(ffvals)];
                    
                    
                    
                    % --- MUSC
                    inds=TvalsMUSC_round==k;
                    ffvals=FFvalsMUSC(inds);
                    cvMUSC=std(ffvals)/abs(mean(ffvals));
                    
                    %                    ffvals_MinusMean=ffvals-mean(ffvals);
                    %                    ffvals_MinusDayMean_AllDays_MUSC=[ffvals_MinusDayMean_AllDays_MUSC ffvals_MinusMean];
                    
                    ffvals_DivideDayMean_AllDays_MUSC=[ffvals_DivideDayMean_AllDays_MUSC ffvals/mean(ffvals)];
                    
                    
                    % ====== save cv vals
                    cvPBS_alldays=[cvPBS_alldays cvPBS];
                    cvMUSC_alldays=[cvMUSC_alldays cvMUSC];
                    
                end
                
                % === get 1) mean of day CVs, and 2) CV over all days
                % (detrended)
                %                 MeanOfDayCVs_PBS=mean(cvPBS_alldays);
                %                 MeanOfDayCVs_MUSC=mean(cvMUSC_alldays);
                
                CVofAllDays_UsingValsDividedByDayMean_PBS=std(ffvals_DivideDayMean_AllDays_PBS);
                CVofAllDays_UsingValsDividedByDayMean_MUSC=std(ffvals_DivideDayMean_AllDays_MUSC);
                
                % ratio of MUSC CV to PBS CV, and whether that is
                % significant
                cvRatio_MUSCoverPBS_usingAllVals=CVofAllDays_UsingValsDividedByDayMean_MUSC/CVofAllDays_UsingValsDividedByDayMean_PBS;
                [~, p]=vartest2(ffvals_DivideDayMean_AllDays_PBS, ffvals_DivideDayMean_AllDays_MUSC, 'tail', 'right');
                if DispEachSylCVpval==1
                    disp([birdname '-' exptname '-' syl ': ' num2str(cvRatio_MUSCoverPBS_usingAllVals), '; p=' num2str(p)]);
                end
                %                 plot(CVofAllDays_UsingValsDividedByDayMean_MUSC, CVofAllDays_UsingValsDividedByDayMean_PBS, 'o');
                
                
                
                % ===== OUTPUT DATA
                % --- cv related stuff
                cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS=[cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS cvRatio_MUSCoverPBS_usingAllVals];
                cvRatio_pvalue_UsingAllVals_ALLSYLS=[cvRatio_pvalue_UsingAllVals_ALLSYLS p];
                
                CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS CVofAllDays_UsingValsDividedByDayMean_PBS];
                CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS CVofAllDays_UsingValsDividedByDayMean_MUSC];
                
                cvPBS_alldays_ALLSYLS=[cvPBS_alldays_ALLSYLS cvPBS_alldays];
                cvMUSC_alldays_ALLSYLS=[cvMUSC_alldays_ALLSYLS  cvMUSC_alldays];
                % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                
                
                % --- other stuff
                Y_PosRelTarg=[Y_PosRelTarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ];
                Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                Y_presimilar=[Y_presimilar SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl];
                %                 Y_AcousticDist=[Y_AcousticDist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore];
                %                 targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
                %                 Y_Corr=[Y_Corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl)];
                
                Yexpt=[Yexpt exptname(end-3:end)];
                Ybird=[Ybird birdname(1:4)];
                
                Expt_count_all=[Expt_count_all expt_count];
                
                
                Y_Generalization_MP=[Y_Generalization_MP Generalization_MP];
                Y_Generalization_AFP=[Y_Generalization_AFP Generalization_AFP];
                Y_Generalization_Learn=[Y_Generalization_Learn Generalization_Learn];
                
                TargLearnDirAll=[TargLearnDirAll targlearndir];
                
            end
            
            
            % ================= Flip sign if learning at targsyl is negative
%             if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
%                 Y_FFmean_pbs=-1.*Y_FFmean_pbs;
%                 Y_FFmean_musc=-1.*Y_FFmean_musc;
%                 Y_AFP_bias=-1.*Y_AFP_bias;
%             end
            
            % ========= Normalize by targsyl if desired (PBS learning
            % by taergsyl)
            if norm_by_targsyl==1
                learning_by_targ=Y_FFmean_pbs(Y_istarg==1)*TargLearnDirAll;
                
                Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
                Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
                Y_AFP_bias=Y_AFP_bias./learning_by_targ;
            end
            
            
            
            % ============================ COLLECT DATA TO PLOT FOR ALL
            % EXPERIMENTS
            if any(~isnan(Y_FFmean_pbs)) % if any are not nan.
                
                % -- cv stuff
                cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS cvRatio_MUSCoverPBS_usingAllVals_ALLSYLS];
                cvRatio_pvalue_UsingAllVals_ALLEXPTS=[cvRatio_pvalue_UsingAllVals_ALLEXPTS cvRatio_pvalue_UsingAllVals_ALLSYLS];
                
                CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_PBS_ALLSYLS];
                CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLSYLS];
                
                cvPBS_alldays_ALLEXPTS=[cvPBS_alldays_ALLEXPTS cvPBS_alldays_ALLSYLS];
                cvMUSC_alldays_ALLEXPTS=[cvMUSC_alldays_ALLEXPTS  cvMUSC_alldays_ALLSYLS];
                
                
                % -- other stuff
                LearningPBS_all=[LearningPBS_all Y_FFmean_pbs];
                MPbias_all=[MPbias_all Y_FFmean_musc];
                AFPbias_all=[AFPbias_all Y_AFP_bias];
                SimDiff_all=[SimDiff_all Y_similar_diff];
                TargStatus_all=[TargStatus_all Y_istarg];
                PreSylSimilar_all=[PreSylSimilar_all Y_presimilar];
                
                Generalization_AFP_all=[Generalization_AFP_all Y_Generalization_AFP];
                Generalization_MP_all=[Generalization_MP_all Y_Generalization_MP];
                Generalization_Learn_all=[Generalization_Learn_all Y_Generalization_Learn];
                
                
                
                Yexpt_all=[Yexpt_all Yexpt];
                Ybird_all=[Ybird_all Ybird];
                Y_PosRelTarg_All=[Y_PosRelTarg_All Y_PosRelTarg];
                
                
                
                expt_count=expt_count+1;
                
            end
        end
    end
end



%% PRE, POST, ETC

lt_figure; hold on;

% ========== [ALL]
lt_subplot(2,3,1); hold on;
xlabel('targ - PRE(1) - PRE(REST) - DIFF MOTIF');
ylabel('AFP bias');
title('line = expt [SAME]');
Yall = [];
Yall_syls = {};
for i=unique(Expt_count_all)
    
   Y = nan(1,4);
   Ycell = cell(1,4);
   
   % == targ
   indsthis = Expt_count_all==i & TargStatus_all==1;
      
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(1) = y;
   Ycell{1} = y;
   
   % === PRE(1)
   indsthis = Expt_count_all==i & TargStatus_all==0 & Y_PosRelTarg_All==-1;
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(2) = mean(y);
   Ycell{2} = y;
   
   % === PRE(REST)
   indsthis = Expt_count_all==i & TargStatus_all==0 & (Y_PosRelTarg_All<-1 & ~isnan(Y_PosRelTarg_All));
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(3) = mean(y);
   Ycell{3} = y;
   
   % === DIFF MOTIF
   indsthis = Expt_count_all==i & TargStatus_all==0 & isnan(Y_PosRelTarg_All);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(4) = mean(y);
   Ycell{4} = y;

%    
%    % === others
%    indsthis = Expt_count_all==i & TargStatus_all==0;
%       
%    x = Y_PosRelTarg_All(indsthis);
%    y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
%       
%    Y(2) = mean(y(~isnan(x))); % same motif.
%    Y(3) = mean(y(isnan(x))); % diff mtoif
   
   plot(1:4, Y, '-ok');
   Yall = [Yall; Y];
   Yall_syls = [Yall_syls; Ycell];
end
lt_plot([1:4]+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
xlim([0 5]);
lt_plot_zeroline;
% [h, p] = ttest(Yall(:,2), Yall(:,3));
% lt_plot_pvalue(p, 'same vs diff motif', 1);


lt_subplot(3,2,2); hold on
title('dat = syl; only expt with all cases (ignore diff)');
tmp = ~cellfun(@isempty, Yall_syls(:,1:3), 'UniformOutput', 1);
% tmp = ~cellfun(@isempty, Yall_syls(:,1:3), 'UniformOutput', 1) | cellfun(@isempty, Yall_syls(:,1:3), 'UniformOutput', 1)
Ytmp = Yall_syls(all(tmp'), :);

for i=1:size(Ytmp,2)
    y = [Ytmp{:,i}];
    plot(i, y, 'ok');
end


%% PRE, POST, REST

lt_figure; hold on;

% ========== [ALL]
lt_subplot(2,3,1); hold on;
xlabel('targ - PRE(1) - REST');
ylabel('AFP bias');
title('line = expt [SAME]');
Yall = [];
Yall_syls = {};
for i=unique(Expt_count_all)
    
   Y = nan(1,4);
   Ycell = cell(1,4);
   
   % == targ
   indsthis = Expt_count_all==i & TargStatus_all==1;
      
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(1) = y;
   Ycell{1} = y;
   
   % === PRE(1)
   indsthis = Expt_count_all==i & TargStatus_all==0 & Y_PosRelTarg_All==-1;
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(2) = mean(y);
   Ycell{2} = y;
   
   % === PRE(REST)
   indsthis = Expt_count_all==i & TargStatus_all==0 & ~(Y_PosRelTarg_All==-1);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(3) = mean(y);
   Ycell{3} = y;

   plot(1:length(Y), Y, '-ok');
   Yall = [Yall; Y];
   Yall_syls = [Yall_syls; Ycell];
end
lt_plot([1:length(Y)]+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
xlim([0 5]);
lt_plot_zeroline;
% [h, p] = ttest(Yall(:,2), Yall(:,3));
% lt_plot_pvalue(p, 'same vs diff motif', 1);
% 
% 
% lt_subplot(3,2,2); hold on
% title('dat = syl; only expt with all cases (ignore diff)');
% tmp = ~cellfun(@isempty, Yall_syls(:,1:3), 'UniformOutput', 1);
% % tmp = ~cellfun(@isempty, Yall_syls(:,1:3), 'UniformOutput', 1) | cellfun(@isempty, Yall_syls(:,1:3), 'UniformOutput', 1)
% Ytmp = Yall_syls(all(tmp'), :);
% 
% for i=1:size(Ytmp,2)
%     y = [Ytmp{:,i}];
%     plot(i, y, 'ok');
% end



%% same, diff, motif

lt_figure; hold on;


% ========== [ALL]
lt_subplot(2,3,1); hold on;
xlabel('targ - same motif - diff motif');
ylabel('AFP bias');
title('line = expt [SAME]');
getsame = 1; % 1 for same, 0 for diff;
Yall =[];
for i=unique(Expt_count_all)
    
   Y = nan(1,3);
   
   % == targ
   indsthis = Expt_count_all==i & TargStatus_all==1;
      
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(1) = y;
   
   % === others
   indsthis = Expt_count_all==i & TargStatus_all==0 & SimDiff_all==getsame;
      
   x = Y_PosRelTarg_All(indsthis);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
      
   Y(2) = mean(y(~isnan(x))); % same motif.
   Y(3) = mean(y(isnan(x))); % diff mtoif
   
   plot(1:3, Y, '-ok');
   
   Yall = [Yall; Y];
end
lt_plot([1:3]+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;
[h, p] = ttest(Yall(:,2), Yall(:,3));
lt_plot_pvalue(p, 'same vs diff motif', 1);


% ========== [ALL]
lt_subplot(2,3,2); hold on;
xlabel('targ - same motif - diff motif');
ylabel('AFP bias');
title('line = expt [DIFF]');

getsame = 0; % 1 for same, 0 for diff;
Yall =[];
for i=unique(Expt_count_all)
    
   Y = nan(1,3);
   
   % == targ
   indsthis = Expt_count_all==i & TargStatus_all==1;
      
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(1) = y;
   
   % === others
   indsthis = Expt_count_all==i & TargStatus_all==0 & SimDiff_all==getsame;
      
   x = Y_PosRelTarg_All(indsthis);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
      
   Y(2) = mean(y(~isnan(x))); % same motif.
   Y(3) = mean(y(isnan(x))); % diff mtoif
   
   plot(1:3, Y, '-ok');
   
   Yall = [Yall; Y];
end
lt_plot([1:3]+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;
[h, p] = ttest(Yall(:,2), Yall(:,3));
lt_plot_pvalue(p, 'same vs diff motif', 1);


% ========== [ALL]
lt_subplot(2,3,3); hold on;
xlabel('targ - same motif - diff motif');
ylabel('AFP bias');
title('line = expt [ALL]');

Yall =[];
for i=unique(Expt_count_all)
    
   Y = nan(1,3);
   
   % == targ
   indsthis = Expt_count_all==i & TargStatus_all==1;
      
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);   
   Y(1) = y;
   
   % === others
   indsthis = Expt_count_all==i & TargStatus_all==0;
      
   x = Y_PosRelTarg_All(indsthis);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
      
   Y(2) = mean(y(~isnan(x))); % same motif.
   Y(3) = mean(y(isnan(x))); % diff mtoif
   
   plot(1:3, Y, '-ok');
   
   Yall = [Yall; Y];
end
lt_plot([1:3]+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;
[h, p] = ttest(Yall(:,2), Yall(:,3));
lt_plot_pvalue(p, 'same vs diff motif', 1);

%%
lt_figure; hold on;

% LearningPBS_all=[];
% 
% MPbias_all=[];
% AFPbias_all=[];
% SimDiff_all=[];
% TargStatus_all=[];
% PreSylSimilar_all=[];
% Expt_count_all=[];
% Yexpt_all={};
% Ybird_all={};
% Y_PosRelTarg_All=[];
% 
% Generalization_MP_all=[];
% Generalization_AFP_all=[];
% Generalization_Learn_all=[];
% 
% cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[];
% cvRatio_pvalue_UsingAllVals_ALLEXPTS=[];
% 
% CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[];
% CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[];
% 
% cvPBS_alldays_ALLEXPTS={};
% cvMUSC_alldays_ALLEXPTS={};
% 
% TargLearnDirAll=[];

Y_PosRelTarg_All(logical(TargStatus_all)) = 0;
Y_PosRelTarg_All(isnan(Y_PosRelTarg_All)) = -9.5;

lt_figure; hold on;

% ========== [ALL]
lt_subplot(2,3,1); hold on;
xlabel('dist from targ');
ylabel('AFP bias');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i;
   
   x = Y_PosRelTarg_All(indsthis);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
    
   X = [X x];
    Y = [Y y];
end
% == overlay means
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;



lt_subplot(2,3,4); hold on;
xlabel('dist from targ');
ylabel('MP bias');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i;
   
   x = Y_PosRelTarg_All(indsthis);
   y = MPbias_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;


% ========== [SAME]
lt_subplot(2,3,2); hold on;
xlabel('dist from targ');
ylabel('AFP bias');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i & SimDiff_all==1 & TargStatus_all==0;
   
   x = Y_PosRelTarg_All(indsthis);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
% == overlay means
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;



lt_subplot(2,3,5); hold on;
xlabel('dist from targ');
ylabel('MP bias');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i & SimDiff_all==1 & TargStatus_all==0;
   
   x = Y_PosRelTarg_All(indsthis);
   y = MPbias_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;





% ========== [DIFF]
lt_subplot(2,3,3); hold on;
xlabel('dist from targ');
ylabel('AFP bias');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i & SimDiff_all==0 & TargStatus_all==0;
   
   x = Y_PosRelTarg_All(indsthis);
   y = AFPbias_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
% == overlay means
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;



lt_subplot(2,3,6); hold on;
xlabel('dist from targ');
ylabel('MP bias');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i & SimDiff_all==0 & TargStatus_all==0;
   
   x = Y_PosRelTarg_All(indsthis);
   y = MPbias_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;







%% [learning, in pbs]

lt_figure; hold on;

% LearningPBS_all=[];
% 
% MPbias_all=[];
% AFPbias_all=[];
% SimDiff_all=[];
% TargStatus_all=[];
% PreSylSimilar_all=[];
% Expt_count_all=[];
% Yexpt_all={};
% Ybird_all={};
% Y_PosRelTarg_All=[];
% 
% Generalization_MP_all=[];
% Generalization_AFP_all=[];
% Generalization_Learn_all=[];
% 
% cvRatio_MUSCoverPBS_usingAllVals_ALLEXPTS=[];
% cvRatio_pvalue_UsingAllVals_ALLEXPTS=[];
% 
% CVofAllDays_UsingValsDividedByDayMean_PBS_ALLEXPTS=[];
% CVofAllDays_UsingValsDividedByDayMean_MUSC_ALLEXPTS=[];
% 
% cvPBS_alldays_ALLEXPTS={};
% cvMUSC_alldays_ALLEXPTS={};
% 
% TargLearnDirAll=[];

Y_PosRelTarg_All(logical(TargStatus_all)) = 0;
Y_PosRelTarg_All(isnan(Y_PosRelTarg_All)) = -9.5;

lt_figure; hold on;

% ========== [ALL]
lt_subplot(2,3,1); hold on;
xlabel('dist from targ');
ylabel('learning');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i;
   
   x = Y_PosRelTarg_All(indsthis);
   y = LearningPBS_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
    
   X = [X x];
    Y = [Y y];
end
% == overlay means
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;



% ========== [SAME]
lt_subplot(2,3,2); hold on;
xlabel('dist from targ');
ylabel('learning');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i & SimDiff_all==1 & TargStatus_all==0;
   
   x = Y_PosRelTarg_All(indsthis);
   y = LearningPBS_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
% == overlay means
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;



% ========== [DIFF]
lt_subplot(2,3,3); hold on;
xlabel('dist from targ');
ylabel('learning');
title('line = expt');
X = [];
Y = [];
for i=unique(Expt_count_all)
   indsthis = Expt_count_all==i & SimDiff_all==0 & TargStatus_all==0;
   
   x = Y_PosRelTarg_All(indsthis);
   y = LearningPBS_all(indsthis).*TargLearnDirAll(indsthis);
   issame = SimDiff_all(indsthis);
   istarg = TargStatus_all(indsthis);
   
   % --- sort
   [~, indstmp] = sort(x);
   x = x(indstmp);
   y = y(indstmp);
   issame = issame(indstmp);
   istarg = istarg(indstmp);
   
   plot(x, y, '-ok');
  
   % -- plot syl type
   plot(x(issame==1 & istarg==0), y(issame==1 & istarg==0), 'ob');
   plot(x(issame==0), y(issame==0), 'or');
   X = [X x];
    Y = [Y y];
end
% == overlay means
[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
lt_plot(unique(X)+0.3, ymean, {'Errors', ysem});
lt_plot_zeroline;








