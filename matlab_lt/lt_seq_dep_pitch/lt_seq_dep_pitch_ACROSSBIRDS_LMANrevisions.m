function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANrevisions(SeqDepPitch_AcrossBirds, PARAMS, norm_by_targsyl, epochfield_input, UseBaselineForCV, DispEachSylCVpval)
%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

if ~exist('norm_by_targsyl','var');
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

epochfield=epochfield_input;

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
            if Y_FFmean_pbs(Y_istarg==1)<0; % negative learning
                Y_FFmean_pbs=-1.*Y_FFmean_pbs;
                Y_FFmean_musc=-1.*Y_FFmean_musc;
                Y_AFP_bias=-1.*Y_AFP_bias;
            end
            
            % ========= Normalize by targsyl if desired (PBS learning
            % by taergsyl)
            if norm_by_targsyl==1;
                learning_by_targ=Y_FFmean_pbs(Y_istarg==1);
                
                Y_FFmean_pbs=Y_FFmean_pbs./learning_by_targ;
                Y_FFmean_musc=Y_FFmean_musc./learning_by_targ;
                Y_AFP_bias=Y_AFP_bias./learning_by_targ;
            end
            
            
            
            % ============================ COLLECT DATA TO PLOT FOR ALL
            % EXPERIMENTS
            if any(~isnan(Y_FFmean_pbs)); % if any are not nan.
                
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

%% ==== display num expts/birds for both targ and nontarg

% TARG
disp(' ==== TARG SYLS');
disp(['numexpts: ' num2str(max(Expt_count_all))]);
disp(['numbirds: ' num2str(length(unique(Ybird_all)))]);

% NONTARG (SAME)
disp(' ==== SAME NONTARG')
inds=find(TargStatus_all==0 & SimDiff_all==1);

exptnumtmp=[];
birdnametmp={};
for i=1:length(inds)
    ind=inds(i);
    
    expt=Expt_count_all(ind);
    bname=Ybird_all{ind};
    
    exptnumtmp=[exptnumtmp expt];
    birdnametmp=[birdnametmp bname];
end

disp(['n: ' num2str(length(inds))])
disp(['numexpts: ' num2str(length(unique(exptnumtmp)))]);
disp(['numbirds: ' num2str(length(unique(birdnametmp)))]);

%% BAR PLOT (FINAL FIGURE FORM)
%% ==== [ALL NONTARG] [same type, diff type] [NO BALLS] [GOOD]
lt_figure; hold on;

Learning_thresh=-10000;
title('All nontarg');

Ylearn_raw={};
YMP_raw={};

% ------ TARGETS
inds=TargStatus_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];


% ------- SIMILAR
inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% -------- DIFF
inds=TargStatus_all==0 & SimDiff_all==0 & LearningPBS_all>Learning_thresh;

learn_raw=LearningPBS_all(inds);
mp_raw=MPbias_all(inds);

Ylearn_raw=[Ylearn_raw learn_raw];
YMP_raw=[YMP_raw mp_raw];



% ================= PLOT (RAW)
X=1:length(Ylearn_raw);

for i=1:length(X);
    xtmp=[i-0.2, i+0.2];
    ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
    plot(xtmp, ytmp', '-k')
end


% ================== PLOT MEAN
Ylearn_mean=cellfun(@mean, Ylearn_raw);
Ymp_mean=cellfun(@mean, YMP_raw);

Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
Ymp_sem=cellfun(@lt_sem, YMP_raw);


lt_plot_bar(X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
hold on;
lt_plot_bar(X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});


% === overlay reversion magnitude
for i=1:length(Ylearn_raw);
    
    ReversionAll=YMP_raw{i}-Ylearn_raw{i};
    revMean=mean(ReversionAll);
    revSEM=lt_sem(ReversionAll);
    lt_plot_text(i-0.1, 1.3*max(YMP_raw{i}), ['ReverMean(sem)= ' num2str(revMean) '(' num2str(revSEM) ')'], 'k')
    lt_plot_text(i-0.1, 1.4*max(YMP_raw{i}), ['PBS= ' num2str(mean(Ylearn_raw{i})) '(' num2str(lt_sem(Ylearn_raw{i})) ')'], 'k')
    lt_plot_text(i-0.1, 1.5*max(YMP_raw{i}), ['MUSC= ' num2str(mean(YMP_raw{i})) '(' num2str(lt_sem(YMP_raw{i})) ')'], 'k')
    
end

% ==== compare reversion magnitudes (raw reversion)
OutputString=[];
for i=1:length(Ylearn_raw)
    Reversion1=YMP_raw{i}-Ylearn_raw{i};
    for ii=i+1:length(Ylearn_raw);
        Reversion2=YMP_raw{ii}-Ylearn_raw{ii};
        
        p=ranksum(Reversion1, Reversion2);
        
        OutputString=[OutputString ['| ' num2str(i) ' vs. ' num2str(ii) ': ' num2str(p)]];
        
    end
end

lt_plot_annotation(1, OutputString, 'm')


% ==== overlay sample size
for i=1:length(Ylearn_raw);
    
    N=numel(Ylearn_raw{i});
    lt_plot_text(i-0.1, 1.2*max(YMP_raw{i}), ['N=' num2str(N)], 'g')
end


% ============ OVERLAY TEXT OF RATIOS (MP/LEARNING)
for i=1:length(Ylearn_raw);
    
    MPoverLearn=mean(YMP_raw{i})/mean(Ylearn_raw{i});
    lt_plot_text(i+0.1, 1.1*max(YMP_raw{i}), [num2str(100*MPoverLearn, '%3.2g') '%'], 'b')
end




for i=1:length(Ylearn_raw);
    
    p = signrank(Ylearn_raw{i}, YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), ['p=' num2str(p)], 'b', 15);
    elseif p<0.005
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
    elseif p<0.05
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
    else
        lt_plot_text(i, max([Ylearn_raw{i} YMP_raw{i}]), ['p=' num2str(p)], 'b', 12);
        
    end

end

% ---------- LEARNING VS 0
for i=1:length(Ylearn_raw);
   
    p = signrank(Ylearn_raw{i});
    
    if p<0.0005
        lt_plot_text(i-0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i-0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i-0.2, 0, '*', 'c', 15);
    end
end


% ---------- MP VS 0
for i=1:length(YMP_raw);
    
    p = signrank(YMP_raw{i});
    
    if p<0.0005
        lt_plot_text(i+0.2, 0, '***', 'c', 15);
    elseif p<0.005
        lt_plot_text(i+0.2, 0, '**', 'c', 15);
    elseif p<0.05
        lt_plot_text(i+0.2, 0, '*', 'c', 15);
    end

end



Xlabels={'Targets','same-type', 'diff-type'};
set(gca, 'XTick', X);
set(gca, 'XTickLabel',Xlabels);
rotateXLabels(gca, 45)

%% == regression analysis (plots, but no regression here)
% --- temp, plotting as in regression analysis
inds = find(TargStatus_all==0 & SimDiff_all==1); % (nontarg context)
lt_figure; hold on;
ylabel('nontarg, ff'); xlabel('targ, ff');
title('blue = pbs; red = musc');

for i=1:length(inds)
    ind = inds(i);
    
    FFpbs_nontarg = LearningPBS_all(ind);
    FFmusc_nontarg = MPbias_all(ind);
    
    exptcount = Expt_count_all(ind);
    indtmp = TargStatus_all==1 & Expt_count_all==exptcount;
    FFpbs_targ = LearningPBS_all(indtmp);
    FFmusc_targ = MPbias_all(indtmp);
    
    plot(FFpbs_targ, FFpbs_nontarg, 'ob');
    plot(FFmusc_targ, FFmusc_nontarg, 'or');
    line([FFpbs_targ FFmusc_targ], [FFpbs_nontarg FFmusc_nontarg], 'Color', 'k');
end 

line(xlim, [0 0], 'Color', 'k');
line([0 0], ylim, 'Color', 'k');
line([0 200], [0 200], 'Color', 'k');


%% ====== POWER ANALYSIS
% Given the effect size for reversion in the target context, how much power
% do we have (given our sample size) to detect reversion in the nontarget
% context

% == METHOD 1 - unpaired ttest 
% figure out effect size
inds = TargStatus_all==0 & SimDiff_all==1; % (nontarg context)
mu0 = mean(LearningPBS_all(inds)); 
sigma0 = std(LearningPBS_all(inds));
Nactual = sum(inds);

% mean reversion as a percent (targ context)
inds=TargStatus_all==1;
mu0tmp = mean(LearningPBS_all(inds));
mu1tmp = mean(MPbias_all(inds));
Reversion_percent = 1 - (mu1tmp/mu0tmp);

% -- expected mu1 if nontarget shows similar reveersion (percent) to targ
mu1 = (1-Reversion_percent)*mu0;

% --- run through diff sample sizes
Nmax = 100; 
PwrAll = [];
for n = 1:Nmax
pwr = sampsizepwr('t', [mu0, sigma0], mu1, [], n);
PwrAll = [PwrAll pwr];
end

lt_figure; hold on;
ylabel('power (to detect reversion in nontarg)'); xlabel('N');
title('unpaired, assuming effect size is reversion(%) in targ context');
plot(1:Nmax, PwrAll, '-');
line([Nactual Nactual], ylim)



% ==== METHOD 2 - for each nontarg, calc reversion based on percent reversion in
% its own targ context, 
% DO THIS WHILE PARAMETRICALLY VARYING REVERSION ( i..e from actual
% revesrion down to a smaller and smaller fraction of actual reveersion)
% NOTE: uses ranksum test (since paired tests will always be significant
% even with low magnitude of reversion, since they are all chagning in same
% direction)
ReversionMultiplierList = 0:0.1:1.5;
PvalueAll = [];
ReversionHzAll = [];
for j=1:length(ReversionMultiplierList)
    
    reversion_multiplier = ReversionMultiplierList(j);
    inds = find(TargStatus_all==0 & SimDiff_all==1); % (nontarg context)
    
    FFmusc_nontarg_expected_ALL = [];
    
    for i=1:length(inds)
        ind = inds(i);
        
        FFpbs_nontarg = LearningPBS_all(ind);
        FFmusc_nontarg = MPbias_all(ind);
        
        exptcount = Expt_count_all(ind);
        indtmp = TargStatus_all==1 & Expt_count_all==exptcount;
        FFpbs_targ = LearningPBS_all(indtmp);
        FFmusc_targ = MPbias_all(indtmp);
        
        % reversion at target for this experiment
        reversion_percent = (FFpbs_targ - FFmusc_targ)/FFpbs_targ;
        reversion_percent = reversion_percent*reversion_multiplier;
        
        % expected nontarg musc, given reversion at targ
        FF_musc_nontarg_expected = (1-reversion_percent)*FFpbs_nontarg;
        
        FFmusc_nontarg_expected_ALL = [FFmusc_nontarg_expected_ALL FF_musc_nontarg_expected];
    end
    
    % ----- If reversion occured at same percent (as in targ) would we detect
    % reversion?
    inds = find(TargStatus_all==0 & SimDiff_all==1); % (nontarg context)
    p = ranksum(LearningPBS_all(inds), FFmusc_nontarg_expected_ALL);
    reversion_hz = mean(LearningPBS_all(inds)) - mean(FFmusc_nontarg_expected_ALL);
    
    PvalueAll = [PvalueAll p];
    ReversionHzAll = [ReversionHzAll reversion_hz];
end

lt_figure; hold on;
lt_subplot(2,1,1); hold on;
xlabel('simulated reversion (multiple of actual reversion in targ');
ylabel('log10(p) from rank sum, nontarg context, simulated');
plot(ReversionMultiplierList, log10(PvalueAll));

lt_subplot(2,1,2); hold on;
xlabel('simulated reversion (multiple of actual reversion in targ');
ylabel('reversion, hz, nontarg context, simulated');
plot(ReversionMultiplierList, ReversionHzAll);


% ===== METHOD 3
% 1) estimate distribtuion of reversion (across all targ cont)

% 2) sample from that distribution, and calculate reversion for each
% nontarget. 
% -- do this by taking fraction to simulate revrsion, then adding on noise
% (in hz) - make the noise distributiont based on the actual standard
% deviation of reversion for nontarget context.

% 3) what is p-value and effect size?

% 4) repeat steps 2 and 3 to get distribution of effect sizes


%% === COMPARE TO WARREN, ANDALMAN AND FEE.

% ==== Andalman
% Fig 3
% plotting deviation from baseline during inactivation vs. during PBS.
% Deviations are mostly between 0 and 80. In almost all cases can see
% reversion. With sample size similar to mine should definitely be able to
% see reversion. Further, their data show individual days, which shoudl
% have more variance. 
% Note: they are pushing up and down so AFP might be more active

% Fig 4
% within day learned pitch and AFP bias. All are between 0 and 50hz
% learning. In almost all cases see reversion. With my sample size would
% definitely see reversion. 

% ===== Warren, 2011
% -- NOTHING - all inactivation days seem to be at poitn of greater
% learning

% ===== Warren, Charlesworth, stim manuscript
% -- within single experiment, see 2 examples for 30-50hz learning and
% reversion.



%% === LOOK AT EARLY DAYS FOR REVERSION IN TARGET WHEN LOWER LEARNING

NumBirds = length(SeqDepPitch_AcrossBirds.birds);
FFpbsAll = [];
FFmuscAll = [];
for i=1:NumBirds
    
    numexpts = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'final_extracted_window')
            disp('SKIPPED');
            continue
        end
        
        % === TARG STATS
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        LastGoodDay = max(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds);
        FirstGoodDay = ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
    
        FFmusc = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow(FirstGoodDay:LastGoodDay);
        FFpbs = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).meanFF_DevFromBase_WithinTimeWindow(FirstGoodDay:LastGoodDay);
        
        dayinds = ~isnan(FFmusc);
        FFmusc = FFmusc(dayinds);
        FFpbs = FFpbs(dayinds);
        
        FFpbsAll = [FFpbsAll FFpbs];
        FFmuscAll = [FFmuscAll FFmusc];
    end
end

lt_figure; hold on;
title('All inactivation days during single dir (all expts)');
lt_plot_45degScatter(FFpbsAll, FFmuscAll, 'k',1);
xlabel('PBS'); ylabel('MUSC');



%% === MATCH TARG/NONTARG, SHOW REVERSION STILL DIFFERENT.

% parametrically vary minimum and maximum learning value - will plot 2d,
% with each "cell" showing one comparison between targ and nontarg

lt_figure; hold on;
LearningMinList = -5:15:115;
LearningMaxList = -5:15:130;
for ll = 1:length(LearningMinList)
    learningMin = LearningMinList(ll);
    
    count = 0;

    for lll = (ll+1):length(LearningMaxList)
        learningMax = LearningMaxList(lll);
        
        count = count+2;
%         lt_subplot(length(LearningMaxList), length(LearningMinList), count); hold on;
        
        Ylearn_raw={};
        YMP_raw={};
        
        % ------ TARGETS
        inds=TargStatus_all==1 & LearningPBS_all>learningMin & LearningPBS_all < learningMax;
        
        learn_raw=LearningPBS_all(inds);
        mp_raw=MPbias_all(inds);
        
        Ylearn_raw=[Ylearn_raw learn_raw];
        YMP_raw=[YMP_raw mp_raw];
        
        
        % ------- SIMILAR
        inds=TargStatus_all==0 & SimDiff_all==1 & LearningPBS_all>learningMin & LearningPBS_all < learningMax;
        
        learn_raw=LearningPBS_all(inds);
        mp_raw=MPbias_all(inds);
        
        Ylearn_raw=[Ylearn_raw learn_raw];
        YMP_raw=[YMP_raw mp_raw];
        
        X=1:length(Ylearn_raw);
        
        
        if length(Ylearn_raw)<2
            continue
        end
    lt_subplot(length(LearningMinList), 1, ll); hold on;
        ylim([-10 150]);
                line([count+0.5 count+0.5], ylim, 'Color','b');

                % === PL;OT
        for i=1:length(X);
            xtmp=[i-0.2, i+0.2];
            ytmp=[Ylearn_raw{i}' YMP_raw{i}'];
            plot(count+xtmp, ytmp', '-k')
        end
        
        
        % ================== PLOT MEAN
        Ylearn_mean=cellfun(@mean, Ylearn_raw);
        Ymp_mean=cellfun(@mean, YMP_raw);
        
        Ylearn_sem=cellfun(@lt_sem, Ylearn_raw);
        Ymp_sem=cellfun(@lt_sem, YMP_raw);
        
        
        lt_plot_bar(count+X-0.2, Ylearn_mean, {'Errors', Ylearn_sem, 'BarWidth', 0.35});
        hold on;
        lt_plot_bar(count+X+0.2, Ymp_mean, {'Errors', Ymp_sem, 'Color', 'r',  'BarWidth', 0.35});
        
        % === REVERSION SIGNIFICNAT>?
        for i=1:length(Ylearn_raw);
            
            p = signrank(Ylearn_raw{i}, YMP_raw{i});
            
            if p<0.0005
                lt_plot_text(count+i, max([Ylearn_raw{i} YMP_raw{i}]), ['p=' num2str(p)], 'b', 15);
            elseif p<0.005
                lt_plot_text(count+i, max([Ylearn_raw{i} YMP_raw{i}]), '**', 'b', 15);
            elseif p<0.05
                lt_plot_text(count+i, max([Ylearn_raw{i} YMP_raw{i}]), '*', 'b', 15);
            elseif p<0.2
                lt_plot_text(count+i, max([Ylearn_raw{i} YMP_raw{i}]), ['p=' num2str(p)], 'b', 12);
            end            
        end
        
        % --- ASK IF REVERSION IS DIFF BETWEEN THE TWO CONTEXTS
                p = ranksum(Ylearn_raw{1}-YMP_raw{1}, Ylearn_raw{2}-YMP_raw{2});
                if p<0.15;
                    lt_plot_text(count+1.5, 1.2*max(Ylearn_raw{1}), ['(rever)p=' num2str(p, '%3.2g')], 'r')
                end
                
                
                % --- ask if the learning is diff between contexts
                 p = ranksum(Ylearn_raw{1}, Ylearn_raw{2});
                if p<0.15;
                    lt_plot_text(count+1.5, 1.3*max(Ylearn_raw{1}), ['(learn)p=' num2str(p, '%3.2g')], 'm')
                end
               
                
                % --- LINE FOR MIN AND MAX
                line([count+0.5 count+2.5], [learningMin learningMin], 'Color', 'm');
                line([count+0.5 count+2.5], [learningMax learningMax], 'Color', 'm');
                
                
    end
    axis tight
end




