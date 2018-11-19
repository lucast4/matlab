function COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct)

%% lt 10/12/18 -
numbirds = length(COHSTRUCT.bird);
for i=1:numbirds
    numexpt = length(COHSTRUCT.bird(i).experiment);
    for ii=1:numexpt
        numsets = length(COHSTRUCT.bird(i).experiment(ii).setnum);
        for ss=1:numsets
            
            nummotifs = length(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif);
            
            % ====== get list of neurons, bregions, and chans for this
            % dataset
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            
            for mm=1:nummotifs
                Chanpairs = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs;
                numchanpairs = size(Chanpairs, 1);

                bregionpairs_sorted = {};
                bregionpairs_unsorted = {};
                for cc=1:numchanpairs
                    chansthis = Chanpairs(cc,:);
                    bregionsthis = {Bregionlist{Chanlist==chansthis(1)}, ...
                        Bregionlist{Chanlist==chansthis(2)}}; assert(length(bregionsthis)==2);
                    
                    bregionpairs_unsorted = [bregionpairs_unsorted; ...
                        [bregionsthis{1} '-' bregionsthis{2}]];
                    
                    % ==== get sorted vbersion
                    bregionsthis = sort(bregionsthis);
                    bregionpairs_sorted = [bregionpairs_sorted; ...
                        [bregionsthis{1} '-' bregionsthis{2}]];
                    
                end
                
                if ~all(strcmp(bregionpairs_sorted, bregionpairs_unsorted)==1)
                    keyboard
                end
                
                
                COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).bregionpairs_sorted = bregionpairs_sorted;
                COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).bregionpairs_unsorted = bregionpairs_sorted;
            end
        end
    end
end

