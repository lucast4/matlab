function COHSTRUCT = lt_neural_Coher_FixCohMatSize(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct, ...
    tlengthdesired, fflengthdesired)

%% lt 10/13/18 - makes sure size of cohmats are identical
% if smaller, then pads with nans. ASSUMES that first timebins are aligned.
% this is reasonable (error will be +/- the size of time bins, since
% coherence is performed on unaligned data, and then they are aligned ...
% if larger, then clips to include only appropriate dat.
numbirds = length(COHSTRUCT.bird);
tbin1All = [];
for i=1:numbirds
    numexpt = length(COHSTRUCT.bird(i).experiment);
    for ii=1:numexpt
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
                
                CohCell = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Coh_ChpairByTrial;
                Chanpairs = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs;
                tbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).t_relons;
                ffbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).ffbins;
                
                if isempty(CohCell)
                    continue
                end
                % ========= go thru all chan pairs and trials and make sure
                % cohmat is correct size
                numchanpairs = size(CohCell,1);
                ntrials = size(CohCell,2);
                tbin1All = [tbin1All; tbins(1)];
                
                for cc=1:numchanpairs
                    
                    % =================== COLLECT COHEROGRAMS ACROSS TRIALS
                    % --------- 1) COLLECT COHEROGRAMS FOR THIS PAIR
                    %                     dim1 = size(CohCell{1},1);
                    %                     dim2 = size(CohCell{1},2);
                    %
                    % --------------- CONVERT CELLS TO A MATRIC
                    % use for loop because soemtime different size..
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
                        if size(CohCell{cc,tt},1)<tlengthdesired | ...
                                size(CohCell{cc, tt}, 2)<fflengthdesired
                            
                            % -- then stcik in, laving edge with nothing
                            tend = size(CohCell{cc,tt},1);
                            ffend = size(CohCell{cc,tt},2);
                            
                            cohmat = nan(tlengthdesired, fflengthdesired);
                            cohmat(1:tend, 1:ffend) = CohCell{cc, tt};
                            % -- stick back in
                            CohCell{cc, tt} = cohmat;
                        elseif size(CohCell{cc,tt},1)>tlengthdesired | ...
                                size(CohCell{cc, tt}, 2)>fflengthdesired
                            
                            cohmat = CohCell{cc, tt}(1:tlengthdesired, 1:fflengthdesired);
                            % -- stick back in
                            CohCell{cc, tt} = cohmat;
                        end
                    end
                end
                % ======== UPDATE STRUCTURE
                COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Coh_ChpairByTrial ...
                    = CohCell;
            end
        end
    end
end