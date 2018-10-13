function COHSTRUCT = lt_neural_Coher_ExtrCohStruct(MOTIFSTATS_pop, SummaryStruct)
%% lt 10/12/18 - extract structure holding coherence for all motifs

motif_predur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
motif_postdur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_postdur;
PrePostRelSameTime = MOTIFSTATS_Compiled.birds(1).Params_regexp.preAndPostDurRelSameTimept;
assert(~isempty(PrePostRelSameTime), 'then I have to check what empty means by default...');

%%

BirdsToPlot = {'pu69wh78', 'wh44wh39'};
% SetsToSkip = {'1-2-2'};
SetsToSkip = {};

COHSTRUCT = struct;
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
            
            
            % ================== FOR THIS SET, COLLECT ALL COHERENCE DATA
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_neurons{ss};
            Fnamelist = MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_songfiles{ss};
            %     dirname = SummaryStruct.birds(i).neurons(Neurlist(1)).dirname;
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            Chanlist_sort = sort(Chanlist);
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            
            % ----- only one chan, then skip
            if length(unique(Chanlist))==1
                continue
            end
            
            % --- get dirname
            dirname_main = {};
            for j=1:length(Neurlist)
                dname = fileparts(SummaryStruct.birds(i).neurons(Neurlist(j)).dirname);
                dirname_main = [dirname_main; dname];
            end
            dirname_main = unique(dirname_main); assert(length(dirname_main)==1);
            dirname_main = dirname_main{1};
            
            % ---- for coherence folder
            dircoh = [dirname_main '/COHERENCE'];
            
            
            % =============== FOR EACH MOTIF, COLLECT ALL TRIALS FOR ALL CHAN PAIRS
            for mm=1:nummotifs
                disp(['motif ' num2str(mm)]);
                if isempty(dat.motif(mm).SegExtr_neurfakeID)
                    continue
                end
                
                % --- segextract is shared across channels
                segextract = dat.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                neurthis = dat.motif(mm).SegExtr_neurfakeID(1).neurID_orig;
                
                
                % =============== EXTRACT COHERENCE FOR THIS MOTIF ACROSS ALL CHANS
                % AND TRIALS
                [CohAllTrials, Chanpairs, t_relons, ffbins] = lt_neural_Coher_GetAllMotifsChans(...
                    SummaryStruct, i, neurthis, segextract, Chanlist, motif_predur, motif_postdur,PrePostRelSameTime);
                
                % ============== COLLECT
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Coh_ChpairByTrial = CohAllTrials;
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).Chanpairs = Chanpairs;
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).t_relons = t_relons;
                COHSTRUCT.bird(i).experiment(ee).setnum(ss).motif(mm).ffbins = ffbins;
            end
        end
    end
end

