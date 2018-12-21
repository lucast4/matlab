function lt_seq_dep_pitch_ACROSSBIRDS_LMANmotif(SeqDepPitch_AcrossBirds, PARAMS)

%% lt 12/19/18 - good code to plto each experiment all syls rel to motif position
%% REMOVE SYLS THAT SHOULD NOT BE ANALYZED (I.E O/L WITH WN, for non-catch analyses)

disp('--');
disp('Removing syllables that should not be analyzed - i.e. WN overlap, since not using catch songs. REMOVED:');
NumBirds = length(SeqDepPitch_AcrossBirds.birds);

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


%% [GOOD] learning vs. mtoif poostion (and musc, pbs)
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments
        
        % ==== ONE FIGURE PER EXPT
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname '-' exptname]);
        
        % ========== COLLECT AND PLOT ONE BY ONE EACH EPOCH (E.G. CONSOLID START...)
        EpochNameList={'final_extracted_window'};
        for k=1:length(EpochNameList)
            epochfield=EpochNameList{k};
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC, epochfield);
                
                
                Y_FFmean_pbs=[];
                Y_FFmean_musc=[];
                Y_FFsem_pbs=[];
                Y_FFsem_musc=[];
                Y_syls={};
                Y_similar_diff=[];
                Y_istarg=[];
                Y_AFP_bias=[];
                
                
                for j=1:length(SylsUnique)
                    syl=SylsUnique{j};
                    
                    % ===== COLLECT DATA - for each syl in order, get learning (PBS and
                    % MUSC)
                    FF_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_pbs; % mean using rends across days
                    FF_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).meanFF_musc;
                    
                    FFsem_PBS=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_pbs;
                    FFsem_MUSC=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).(syl).semFF_musc;
                    
                    
                    % ===== OUTPUT DATA
                    Y_FFmean_pbs=[Y_FFmean_pbs FF_PBS];
                    Y_FFmean_musc=[Y_FFmean_musc FF_MUSC];
                    Y_FFsem_pbs=[Y_FFsem_pbs FFsem_PBS];
                    Y_FFsem_musc=[Y_FFsem_musc FFsem_MUSC];
                    Y_syls=[Y_syls, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl];
                    Y_similar_diff=[Y_similar_diff SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ];
                    Y_istarg=[Y_istarg SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target];
                    Y_AFP_bias=[Y_AFP_bias FF_PBS-FF_MUSC];
                end
                
                % ====== PLOT
                %                 lt_subplot(2,2,k); hold on;
                
                % #############################################
                X = 1:length(Y_AFP_bias);
                
                lt_plot_bar(X-0.3, Y_FFmean_pbs, {'Color', 'k', 'Errors', Y_FFsem_pbs, 'BarWidth', 0.25});
                lt_plot_bar(X, Y_FFmean_musc, {'Color', 'r', 'Errors', Y_FFsem_musc, 'BarWidth', 0.25});
                lt_plot_bar(X+0.3, Y_AFP_bias, {'Color', 'w', 'BarWidth', 0.25});
                
                % ===== NOTE TARGET
                YLIM = ylim;
                plot(find(Y_istarg==1), 0.9*(YLIM(2)), 'rd');
                
                indstmp = find(Y_istarg==0 & Y_similar_diff==1);
                if ~isempty(indstmp)
                    plot(indstmp, 0.9*(YLIM(2)), 'bd');
                end
                % ++++++++++++++++++++++++++++++++++++++++++++++ GLOBAL
                set(gca, 'XTick', X);
                set(gca, 'XTickLabel', Y_syls)
                
                %                 % ====== what were the days for this epoch?
                %                 days_used=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epochfield).days_list;
                %
                %                 title([epochfield(6:10) epochfield(end-5:end) '-d' num2str(days_used)]);
                %
                
            end
        end
        
        % ===== GLOBAL
%         lt_subtitle([birdname '-' exptname]);
    end
end

