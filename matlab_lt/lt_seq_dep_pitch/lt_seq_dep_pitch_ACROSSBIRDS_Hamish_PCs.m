function lt_seq_dep_pitch_ACROSSBIRDS_Hamish_PCs(SeqDepPitch_AcrossBirds, PARAMS, ...
        targonly, baseonly)

%% lt 7/27/18 -

% max rends to plot for a given day
Nmax = 30;

%%
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

%%

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts
        
        % =============================== GET SET OF SYLS TO COLLECT
        SylsToCollect = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        DatThis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning;
        
            
            
        
        if isfield(DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'final_extracted_window')
            
        % -- max baseline day
        WNday1 = DatThis.Params.PlotLearning.WNTimeOnInd;
        if baseonly==1
            maxday = WNday1-1;
        elseif baseonly ==0
            % - then plot all days
            maxday = DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds(end);
        end
        
        % ============================= timebins same for all syls
            tbinsize = DatThis.Params.DayRawDat.pc_T{1}(2) - DatThis.Params.DayRawDat.pc_T{1}(1);
            t1 = DatThis.Params.DayRawDat.pc_T{1}(1);
            
            for j=1:length(DatThis.Params.DayRawDat.pc_T)
                assert(DatThis.Params.DayRawDat.pc_T{j}(1) == t1);
                % --- since I am assuming that for all syls t1 is the same
                % ...
            end
                
            hsplots = [];
            for ss = 1:length(SylsToCollect)
                sylthis = SylsToCollect{ss};
                
                istarg = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(sylthis).is_target;
                if targonly==1
                    if istarg==0
                        continue
                    end
                end
                
                % ----------
                birdname = SeqDepPitch_AcrossBirds.birds{i}.birdname;
                exptname = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
                
                fdirthis = ['/bluejay5/lucas/analyses/seqdeppitch/GetPitchCont/' birdname '_' exptname];
                
                
                % ====================================== PLOT
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(sylthis);
                ylabel([birdname '-' exptname]);
                
                for nn=1:maxday
                    daythis = nn;
                    
                    % ======== PBS;
                    PBSorMUSC = 'PBS';
                    
                    % ----------- 1) load dat (PCmat
                    fnamethis = [fdirthis '/syl' sylthis '_' PBSorMUSC '_day' num2str(daythis) '.mat'];
                    if exist(fnamethis, 'file')
                        PCmat = load(fnamethis);
                        PCmat = PCmat.PCmat_tosave;
                      
                        if size(PCmat,1)>Nmax
                            indstmp = randperm(size(PCmat,1), Nmax);
                            PCmat = PCmat(indstmp, :);
                        end
                        
                        T = t1+(0:size(PCmat,2)-1)*tbinsize;
                        plot(T, PCmat, '-', 'Color', [0.7 0.7 0.7]);
                    end
                    
                    
                    % ======== PBS;
                    PBSorMUSC = 'MUSC';
                    
                    % ----------- 1) load dat (PCmat
                    fnamethis = [fdirthis '/syl' sylthis '_' PBSorMUSC '_day' num2str(daythis) '.mat'];
                    if exist(fnamethis, 'file')
                        PCmat = load(fnamethis);
                        PCmat = PCmat.PCmat_tosave_MUSC;
                      
                        if size(PCmat,1)>Nmax
                            indstmp = randperm(size(PCmat,1), Nmax);
                            PCmat = PCmat(indstmp, :);
                        end
                        
                        T = t1+(0:size(PCmat,2)-1)*tbinsize;
                        plot(T, PCmat, '-', 'Color', [0.7 0.2 0.2]);
                    end

                    % -----------------------
                    axis tight;
                end
                
            end
            
            linkaxes(hsplots, 'x');
        end
    end
end
