function LFPSTRUCT = lt_neural_LFP_ExtractStruct(MOTIFSTATS_pop, SummaryStruct, ...
    MOTIFSTATS_Compiled, skipifOnlyOneChan, BirdsToPlot, SetsToSkip)
%% lt 10/29/18 - analogous to lt_neural_Coher_ExtrCohStruct except extracts LFP and them get spectrogram.

% BirdsToPlot = {'pu69wh78', 'wh44wh39'};
% % SetsToSkip = {'1-2-2'};
% SetsToSkip = {};

%%
motif_predur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
extrapad = 0.05; % collect extra (sec) on edges.

%%
LFPSTRUCT = struct;
% ==============
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    disp(birdname);
    
    if ~any(strcmp(birdname, BirdsToPlot))
        continue
    end
    
    numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);
    for ee = 1:numexpts
        disp(['expt ' num2str(ee)]);
        
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_neurons);
        
        % =============== FOR THIS BIRD, COLLECT ALL COHERENCE
        for ss = 1:numsets
            disp(['neuron set ' num2str(ss)]);
            
            if any(strcmp([num2str(i) '-' num2str(ee) '-' num2str(ss)], SetsToSkip))
                % then skip this one
                disp(['Skipping (bird, exp, set): ' [num2str(i) '-' num2str(ee) '-' num2str(ss)]]);
                continue
            end
            
            dat = MOTIFSTATS_pop.birds(i).exptnum(ee).DAT.setnum(ss);
            nummotifs = length(dat.motif);
            
            % ================== FOR THIS SET, COLLECT ALL LFP DATA
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_neurons{ss};
            %             Fnamelist = MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_songfiles{ss};
            %             %     dirname = SummaryStruct.birds(i).neurons(Neurlist(1)).dirname;
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            %             Chanlist_sort = sort(Chanlist);
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            [~, indsort] = sort(Chanlist);
            Chanlist = Chanlist(indsort);
            Bregionlist = Bregionlist(indsort);
            
            
            
            % ----- only one chan, then skip
            if skipifOnlyOneChan==1
                if length(unique(Chanlist))==1
                    continue
                end
            end
            
            % --- get dirname
            dirname_main = {};
            for j=1:length(Neurlist)
                dname = fileparts(SummaryStruct.birds(i).neurons(Neurlist(j)).dirname);
                dirname_main = [dirname_main; dname];
            end
            dirname_main = unique(dirname_main); assert(length(dirname_main)==1);
            dirname_main = dirname_main{1};
            
            % =============== FOR EACH MOTIF, COLLECT ALL TRIALS FOR ALL
            % CHANNELS
            for mm=1:nummotifs
                disp(['motif ' num2str(mm)]);
                if isempty(dat.motif(mm).SegExtr_neurfakeID)
                    continue
                end
                
                % --- segextract is shared across channels
                segextract = dat.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                %                 neurthis = dat.motif(mm).SegExtr_neurfakeID(1).neurID_orig;
                
                
                % ================ COLLECT LFP FOR ALL CHANNELS
                [LFPall, Tall, chanlist2] = lt_neural_QUICK_Segextr_GetLFP(segextract, dirname_main, Chanlist, ...
                    motif_predur, extrapad);
                if size(Chanlist,2)>1
                    assert(all(chanlist2==Chanlist'));
                else
                    assert(all(chanlist2==Chanlist));
                end
                
                % ============== COLLECT
                LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).LFP_chanbytrial = LFPall;
                LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Chanlist = Chanlist;
                LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Bregionlist = Bregionlist;
                LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).t_relons = Tall{1};
                %                 LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Coh_ChpairByTrial = CohAllTrials;
                %                 LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Chanpairs = Chanpairs;
                %                 LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).t_relons = t_relons;
                %                 LFPSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).ffbins = ffbins;
            end
        end
    end
end