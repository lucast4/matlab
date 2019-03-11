close all;
savebase = '/bluejay0/bluejay2/lucas/analyses/neural/LFP/FIGS_PlotEgRaw_v2/05Feb2019_2142';

for i=1:length(SwitchCohStruct.bird)
% for i=3
    for ii=1:length(SwitchCohStruct.bird(i).exptnum)
        for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
            
            %% extract
            
            birdplot = SwitchStruct.bird(i).birdname;
            exptplot = SwitchStruct.bird(i).exptnum(ii).exptname;
            swplot = ss;
            motifplot = []; % [string] leave blank for target
            
            
            % #####################################################
            %             i = find(strcmp({SummaryStruct.birds.birdname}, birdplot));
            %             ii = find(strcmp({SummaryStruct.birds(i).exptnum_pop.exptname}, exptplot));
            %             ss = swplot;
            
            indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
                OUTSTRUCT.istarg==1);
            
            if ~any(indsthis)
                continue
            end
            disp([i ii ss]);
            
            % --- target syl motif number?
            motifnumlist = OUTSTRUCT.motifnum(indsthis);
            motifplotlist = OUTSTRUCT.motifname(indsthis);
            
            [~, indstmp] = unique(motifnumlist);
            
            motifnumlist = motifnumlist(indstmp);
            motifplotlist = motifplotlist(indstmp);
            
            %             if length(mm)>1
            %                 disp(mm);
            %                 disp(motifplot);
            %                 disp('MULTIPLE TARGS, taking first one');
            %                 mm = mm(1);
            %                 indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
            %                     OUTSTRUCT.istarg==1 & OUTSTRUCT.motifnum==mm);
            %                 motifplot = OUTSTRUCT.motifname{indsthis(1)};
            %             end
            %             assert(length(mm)==1, 'multipel targets?');
            %             disp(motifplot)
            
            for j=1:length(motifplotlist)
                mm = motifnumlist(j);
                motifplot = motifplotlist{j};
                
                disp(mm);
                disp(motifplot);
                
                if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum)
                    continue
                end
                
                % ======== 1) EXTRACT RAW DATA (NOT JUST LFP)
                extrapad = 0.05; % seconds, pre and post...
                [DatAll, t_onoff, fs, bregionlist, chanlist_toget] = ...
                    lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, ...
                    SummaryStruct, SwitchCohStruct, birdplot, exptplot, swplot, motifplot, ...
                    extrapad);
                
                
                % ========= initiate save dir for this bird
                sdir_this = [savebase '/' birdplot '_' exptplot '_sw' num2str(swplot) '_' motifplot];
                
                %% plot heat mat (in trial order)
                filt_low = 20;
                filt_hi = 35;
                noiseband = [3100 3400];
                plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
                %                 trialOrder = 'coh';
                trialOrder = 'trials'; % in order of actual trials
                cohwind_f = PARAMS.cohscal_fwind;
                
                lt_neural_LFP_PlotRawHeatmap(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                    i, ii, swplot, mm, filt_low, filt_hi, SwitchCohStruct, ...
                    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
                    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered, trialOrder, ...
                    noiseband, cohwind_f);
                
                mkdir([sdir_this '/alltrials_inorder']);
                cd([sdir_this '/alltrials_inorder']);
                lt_save_all_figs;
                close all;
                
                %% plot heat mat (coherence order)
                filt_low = 20;
                filt_hi = 35;
                noiseband = [3100 3400];
                plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
                trialOrder = 'coh';
                % trialOrder = 'trials'; % in order of actual trials
                cohwind_f = PARAMS.cohscal_fwind;
                
                lt_neural_LFP_PlotRawHeatmap(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                    i, ii, swplot, mm, filt_low, filt_hi, SwitchCohStruct, ...
                    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
                    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered, trialOrder, ...
                    noiseband, cohwind_f);
                
                mkdir([sdir_this '/alltrials_increaseingcoh']);
                cd([sdir_this '/alltrials_increaseingcoh']);
                lt_save_all_figs;
                close all;
                
                
                %% plot some trials. [hi]
                
                % ================== FIND TRIAL TO PLOT
                % trialmode = 'random';
                trialmode = 'highcoh';
                % trialmode = 'lowcoh';
                N = 10;
                
%                 ntrials = size(DatAll{1},2);
                
                % ---- list of all possible trials
                indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
                    OUTSTRUCT.motifnum==mm);
                triallist = [find(OUTSTRUCT.indsbase{indsthis(1)}) find(OUTSTRUCT.indsWN{indsthis(1)})];
                
                % ---- what is vector of coherence
                cohscal = OUTSTRUCT.cohscal(indsthis);
                cohscal = mean(cell2mat(cohscal),1); % take average over all chan pairs
                
                [~, indsort] = sort(cohscal);
                indsort = indsort(ismember(indsort, triallist));
                indslow = indsort(1:N);
                indshi = indsort(end-N+1:end);
                
                
                if strcmp(trialmode, 'random')
                    trialstoget = triallist(randperm(length(triallist), N));
                elseif strcmp(trialmode, 'highcoh')
                    trialstoget = indshi;
                elseif strcmp(trialmode, 'lowcoh')
                    trialstoget = indslow;
                end
                
                % ======================================= PLOT ALL THINGS COMBINED
                % _------------ MAKE SURE TO FIRST HAVE COHEROGRAMS...
                close all;
                filt_low = 20;
                filt_hi = 35;
                % filt_low = 3200;
                % filt_hi = 3300;
                % trialtoplot = 100;
                plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
                lt_neural_LFP_PlotEgRaw_v2(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                    i, ii, swplot, mm, trialstoget, filt_low, filt_hi, SwitchCohStruct, ...
                    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
                    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered);
                
                
                
                
                mkdir([sdir_this '/randomtrials_hiCoh']);
                cd([sdir_this '/randomtrials_hiCoh']);
                lt_save_all_figs;
                close all;
                %% plot some trials. [lo]
                
                % ================== FIND TRIAL TO PLOT
                % trialmode = 'random';
                % trialmode = 'highcoh';
                trialmode = 'lowcoh';
                N = 10;
                
                ntrials = size(DatAll{1},2);
                
                % ---- list of all possible trials
                indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
                    OUTSTRUCT.motifnum==mm);
                triallist = [find(OUTSTRUCT.indsbase{indsthis(1)}) find(OUTSTRUCT.indsWN{indsthis(1)})];
                
                % ---- what is vector of coherence
                cohscal = OUTSTRUCT.cohscal(indsthis);
                cohscal = mean(cell2mat(cohscal),1); % take average over all chan pairs
                
                [~, indsort] = sort(cohscal);
                indsort = indsort(ismember(indsort, triallist));
                indslow = indsort(1:N);
                indshi = indsort(end-N+1:end);
                
                
                if strcmp(trialmode, 'random')
                    trialstoget = triallist(randperm(length(triallist), N));
                elseif strcmp(trialmode, 'highcoh')
                    trialstoget = indshi;
                elseif strcmp(trialmode, 'lowcoh')
                    trialstoget = indslow;
                end
                
                % ======================================= PLOT ALL THINGS COMBINED
                % _------------ MAKE SURE TO FIRST HAVE COHEROGRAMS...
                close all;
                filt_low = 20;
                filt_hi = 35;
                % filt_low = 3200;
                % filt_hi = 3300;
                % trialtoplot = 100;
                plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
                lt_neural_LFP_PlotEgRaw_v2(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                    i, ii, swplot, mm, trialstoget, filt_low, filt_hi, SwitchCohStruct, ...
                    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
                    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered);
                
                
                
                
                mkdir([sdir_this '/randomtrials_LoCoh']);
                cd([sdir_this '/randomtrials_LoCoh']);
                lt_save_all_figs;
                close all;
                %% plot some trials. [any]
                
                % ================== FIND TRIAL TO PLOT
                trialmode = 'random';
                % trialmode = 'highcoh';
                % trialmode = 'lowcoh';
                N = 10;
                
                ntrials = size(DatAll{1},2);
                
                % ---- list of all possible trials
                indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
                    OUTSTRUCT.motifnum==mm);
                triallist = [find(OUTSTRUCT.indsbase{indsthis(1)}) find(OUTSTRUCT.indsWN{indsthis(1)})];
                
                % ---- what is vector of coherence
                cohscal = OUTSTRUCT.cohscal(indsthis);
                cohscal = mean(cell2mat(cohscal),1); % take average over all chan pairs
                
                [~, indsort] = sort(cohscal);
                indsort = indsort(ismember(indsort, triallist));
                indslow = indsort(1:N);
                indshi = indsort(end-N+1:end);
                
                
                if strcmp(trialmode, 'random')
                    trialstoget = triallist(randperm(length(triallist), N));
                elseif strcmp(trialmode, 'highcoh')
                    trialstoget = indshi;
                elseif strcmp(trialmode, 'lowcoh')
                    trialstoget = indslow;
                end
                
                % ======================================= PLOT ALL THINGS COMBINED
                % _------------ MAKE SURE TO FIRST HAVE COHEROGRAMS...
                close all;
                filt_low = 20;
                filt_hi = 35;
                % filt_low = 3200;
                % filt_hi = 3300;
                % trialtoplot = 100;
                plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
                lt_neural_LFP_PlotEgRaw_v2(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                    i, ii, swplot, mm, trialstoget, filt_low, filt_hi, SwitchCohStruct, ...
                    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
                    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered);
                
                
                
                
                mkdir([sdir_this '/randomtrials']);
                cd([sdir_this '/randomtrials']);
                lt_save_all_figs;
                close all;
                
            end
        end
    end
end
