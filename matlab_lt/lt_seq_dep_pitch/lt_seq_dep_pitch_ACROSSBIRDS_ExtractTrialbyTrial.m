function [TrialStruct, Params] = lt_seq_dep_pitch_ACROSSBIRDS_ExtractTrialbyTrial(SeqDepPitch_AcrossBirds, ...
    OnlyExptsWithNoStartDelay, DayWindow, onlyIfSignLearn, useHandLabSame)
%% 7/31/18 - lt, fix
% fixing slight issue where for LMAN experiments collected deviation FF,
% while non-LMAN collected raw FF. did not affect any subsequent analyses,
% because they all used deviation from baseline.

usefix =1; % leave at 1 - default. sets all to use raw FF

%%

% OnlyExptsWithNoStartDelay= 0 ;
% TakeIntoAccountStartDelay = 1;
% DayWindow = [-2 4]; % [-2 4] mean 2 base days and 1st 4 learning days

TakeIntoAccountStartDelay = 1;
%%


Params.OnlyExptsWithNoStartDelay= OnlyExptsWithNoStartDelay;
Params.TakeIntoAccountStartDelay = TakeIntoAccountStartDelay;
Params.DayWindow = DayWindow; %


%% =========== only significnat learning?
if onlyIfSignLearn==1
    filter='learning_metric';
    [SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(...
        SeqDepPitch_AcrossBirds, filter);
    % 40.34hz
    % 0.8z
end

%% lt 8/18/17 - perform cross correlation analysis, to look at lag of generalization

NumBirds = length(SeqDepPitch_AcrossBirds.birds);

for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    TrialStruct.birds(i).birdname = birdname;
    
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        TrialStruct.birds(i).exptnum(ii).exptname = exptname;
        
        % ===== SKIP IF HAD DELAY FROM WN TO EXPT START
        if OnlyExptsWithNoStartDelay==1
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode>0
                disp(['SKIPPED (delay after WN start): ' birdname '-' exptname]);
                TrialStruct.birds(i).exptnum(ii).sylnum = [];
                continue
            end
        end
        
        
        % ------ what data to use?
        LMANinactivated = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated;
        if LMANinactivated==1
            if usefix==1
                datafield_ff='FFvals_WithinTimeWindow';
            elseif usefix==0
                datafield_ff='FFvals_DevFromBase_WithinTimeWindow';
            end
            datafield_tvals='Tvals_WithinTimeWindow';
        else
            datafield_ff='FFvals';
            datafield_tvals='Tvals';
        end
        
        %% ----------------- what were good learning days?
        % ==================== method 1
        numEmptyDays = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
        LearningDaysUsed = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        
        assert(all(LearningDaysUsed == ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.dayIndsUsed), 'asfasdf');
        assert(WNday1+numEmptyDays + 2 == LearningDaysUsed(1), 'asfsdaf');
        
        % ==== for each laerning day go to raw dat folder and extract num
        % renditions
        LearningDays = [WNday1+numEmptyDays:LearningDaysUsed(2)];
        assert(length(LearningDays) == 4, 'afasdf');
        
        % ==================== method 2
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        
        % -- baseline days
        base1=WNday1+DayWindow(1);
        baseEnd=WNday1-1;
        
        % -- WN days
        if TakeIntoAccountStartDelay==1
            daysDelay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
            if daysDelay>0
                disp(['ADDED ' num2str(daysDelay) ' days for: ' birdname '-' exptname]);
            end
        else
            daysDelay=0;
        end
        
        WN1=WNday1+daysDelay;
        WNend=WN1+DayWindow(2)-1;
        
        DaysToExtract=[base1:baseEnd WN1:WNend];
        BaseDays=[base1:baseEnd];
        WNDays=[WN1:WNend];
        
        % ------
        targlearndir = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        % ================ COLLECT SOME GENERAL INFOR FOR THIS EXPT
        assert(~isempty(BaseDays),'asdfasd');
        
        TrialStruct.birds(i).exptnum(ii).BaseDays = BaseDays;
        TrialStruct.birds(i).exptnum(ii).WNDays = WNDays;
        TrialStruct.birds(i).exptnum(ii).WNday1 = WNday1;
        TrialStruct.birds(i).exptnum(ii).SylsUnique = SylsUnique;
        TrialStruct.birds(i).exptnum(ii).numEmptyDays = numEmptyDays;
        TrialStruct.birds(i).exptnum(ii).LMANinactivated =LMANinactivated;
        TrialStruct.birds(i).exptnum(ii).targlearndir =targlearndir;
        TrialStruct.birds(i).exptnum(ii).targsyl =targsyl;
        
        
        %% Collect trial by trail pitch and times for each day and syl
        
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            
            istarget=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
            if useHandLabSame ==0
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            elseif useHandLabSame==1
                similar = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab;
            end
            
            
            % --- COLLECT
            FFvals = [];
            Tvals = [];
            for jj=DaysToExtract
                ffvals = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_ff){jj};
                tvals = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).(datafield_tvals){jj};
                
                if LMANinactivated ==0
                    ffvals = cell2mat(ffvals);
                    tvals = cell2mat(tvals);
                end
                
                % ---- COLLECT
                FFvals = [FFvals; ffvals'];
                Tvals = [Tvals; tvals'];
            end
            
            % --- convert Tvals to days
            firstday = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
            tmp =   lt_convert_EventTimes_to_RelTimes(firstday, Tvals);
            Tvals_datenum = Tvals;
            Tvals = tmp.FinalValue;
            
            missingsomedat = 0;
            if length(DaysToExtract) > length(unique(floor(Tvals))')
                missingsomedat = 1;
            end
            
            % ------------ SORT BY TVALS
            [~, inds] = sort(Tvals);
            Tvals = Tvals(inds);
            FFvals = FFvals(inds);
            Tvals_datenum = Tvals_datenum(inds);
            
            % ---- OUT
            TrialStruct.birds(i).exptnum(ii).sylnum(j).syl = syl;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals = Tvals;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals_datenum = Tvals_datenum;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).FFvals = FFvals;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_istarget = istarget;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_similar = similar;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_SylDimensions = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl);
            TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_missingsomedat = missingsomedat;
            
            % ==== sanity check, make sure this is same value that
            % extracted learning value
            
            if DayWindow(1)==-2 & DayWindow(2)==4
                % last 2 days
                ffmeanWN = mean(FFvals(Tvals>=WNDays(end-1)));
                ffmeanBase = mean(FFvals(Tvals< BaseDays(2)+1));
                
                learnval = (ffmeanWN - ffmeanBase);
                
                assert(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean - learnval < 10, 'assfasd');
            end
            
        end
    end
end

