function [All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, ...
    All_chanpair, All_bregionpair, All_bregionpair_alphaorder, ...
    All_tbins, All_ffbins] = lt_neural_Coher_Summ_Extr(COHSTRUCT, MOTIFSTATS_pop, ...
    SummaryStruct, tlengthdesired, fflengthdesired, PARAMS)
%% === lt 10/12/18 - Series of summary analyses - here just does extraction

savebase = ['/bluejay5/lucas/analyses/neural/LFP/PROCESSED/' PARAMS.savemarker];
%%
All_CohgramMean = {};
All_birdnum = [];
All_enum = [];
All_setnum = [];
All_motifname = {};
All_chanpair = [];
All_bregionpair = {};
All_bregionpair_alphaorder = {};

All_tbins = {};
All_ffbins = {};
numbirds = length(COHSTRUCT.bird);
tbin1All = [];
for i=1:numbirds
    numexpt = length(COHSTRUCT.bird(i).experiment);
    birdname = SummaryStruct.birds(i).birdname;
    
    for ii=1:numexpt
        exptid = MOTIFSTATS_pop.birds(i).exptnum(ii).exptname;
        numsets = length(COHSTRUCT.bird(i).experiment(ii).setnum);
        
        disp([num2str(i) '-' num2str(ii)]);
        
        for ss=1:numsets
            
        
            nummotifs = length(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif);
            
            % ====== get list of neurons, bregions, and chans for this
            % dataset
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            
            for mm=1:nummotifs
                motifname = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss).motif(mm).regexpstr;
                
                if isempty(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs)
                    continue
                end
                   
                
                if isfield(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm), 'Coh_ChpairByTrial')
                    CohCell = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Coh_ChpairByTrial;
                else
                    % then load it
                    filename = ...
                        [savebase '/Coh_bird' num2str(i) '_expt' num2str(ii) '_set' num2str(ss) '_mot' num2str(mm) '.mat'];
                    pairstoget = [];
                    CohCell = lt_neural_LFP_loadProcessDat(filename, pairstoget, 1);
                end
                    
                Chanpairs = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs;
                tbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).t_relons;
                ffbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).ffbins;
                
                if isempty(CohCell)
                    continue
                end
                
                tbin1All = [tbin1All; tbins(1)];
                numchanpairs = size(Chanpairs, 1);
                
                for cc=1:numchanpairs
                    chansthis = Chanpairs(cc,:);
                    bregionsthis = {Bregionlist{Chanlist==chansthis(1)}, ...
                        Bregionlist{Chanlist==chansthis(2)}}; assert(length(bregionsthis)==2);
                    
                    % =================== COLLECT COHEROGRAMS ACROSS TRIALS
                    % --------- 1) COLLECT COHEROGRAMS FOR THIS PAIR
                    %                     dim1 = size(CohCell{1},1);
                    %                     dim2 = size(CohCell{1},2);
                    %
                    % --------------- CONVERT CELLS TO A MATRIC
                    % use for loop because soemtime different size..
                    ntrials = size(CohCell,2);
                    cohmat = nan(length(tbins), length(ffbins), ntrials);

                    for tt=1:ntrials
                        % =================== NOTE, HERE IS AD HOC, ASSUMES
                        % THAT T(1) IS THE SAME NO MATTER THE SIZE OF THE
                        % ARRAY. PROBLEM IS THAT i DID NOT SAVE ALL TBINS
                        % FOR ALL TRIALS. BUT AI THINK THAT IS THE CASE
                        % BASED ONT HE CASES IN HWIHC ALL TRIALS HAVE A
                        % DIFFERENT TBIN FROM THE OTHER MOTIFS, AND EVEN
                        % THEN THE T(1) IST HE SAME AS THE OTHER MTOIFS ...
                        % TO DO: FIT ODD-SIZED DATA TO CORRECT POSOTION
                        % EARLY ON IN ANALYSES -I.E. IN EXTACTION OF
                        % COHSTRUCT./
                        if size(CohCell{cc,tt},1)<tlengthdesired
                            % -- then stcik in, laving edge with nothing
                            tend = size(CohCell{cc,tt},1);
                        else
                            tend = tlengthdesired;
                        end
                        cohmat(1:tend, :, tt) = CohCell{cc, tt}(1:tend, 1:fflengthdesired);
%                         if all(size(cohmat(:,:,tt))==size(CohCell{cc,tt}))
%                             cohmat(:,:,tt) = CohCell{cc,tt};
%                         else
%                             disp('skipped! wrong size');
%                         end
                    end
                    
                    % =========== get mean coherogram
                    cohmean = nanmean(cohmat, 3);
                    
                    %                     cohmat = reshape(cell2mat(CohCell(cc, :)), dim1, dim2, []);
                    
                    % ====================
                    tbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).t_relons;
                    ffbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).ffbins;
                    
                    All_tbins = [All_tbins; tbins];
                    All_ffbins = [All_ffbins; ffbins];

                    % ========================== SAVE OUTPUT
                    All_CohgramMean = [All_CohgramMean; cohmean];
                    All_birdnum = [All_birdnum; i];
                    All_enum = [All_enum ; ii];
                    All_setnum = [All_setnum ; ss];
                    All_motifname = [All_motifname; motifname];
                    All_chanpair = [All_chanpair; chansthis];
                    All_bregionpair = [All_bregionpair; bregionsthis];
                    
                    bregionsthis_sort = sort(bregionsthis);
                    All_bregionpair_alphaorder = [All_bregionpair_alphaorder; ...
                        [bregionsthis_sort{1} '-' bregionsthis_sort{2}]];
                end
            end
        end
    end
end

assert(all(diff(tbin1All)<0.001), 'problem with assumption that all are aligned at t1/...');