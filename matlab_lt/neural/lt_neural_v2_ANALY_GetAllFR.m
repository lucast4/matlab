%% lt 6/30/17 - across all birds, expts, neurons, extract raw FR

function [MOTIFSTATS_Compiled] = lt_neural_v2_ANALY_GetAllFR(MOTIFSTATS_Compiled, ...
    RemoveTrialsZeroFR, premotorWind)

%%
NumBirds = length(MOTIFSTATS_Compiled.birds);

%% GET MEAN FR - go thru all experiments

for i=1:NumBirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    
    for ii=1:numexpts
        
        MotifStats =  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        motiflist = MotifStats.params.motif_regexpr_str;
        
        
        % ================ EACH bird, expt, syl
        nummotifs = length(motiflist);
        
        for j=1:nummotifs
            
            disp(['bird' num2str(i) ', exot' num2str(ii) ', motif' num2str(j)]);
            disp('  ');
            % 1) === PLOT LEARNING, OVERLAY OVER TARGET LEARNING. INCLUDE TIME
            % INTERVALS FOR NEURONS
            
            % 2) === FOR EACH NEURON, PLOT CHANGE OVER TIME
            numneurons = length(MotifStats.neurons);
            
            for nn=1:numneurons
                
                
                % ==== For each rendition, get smoothed firing rate
                clustnum = MotifStats.neurons(nn).clustnum;
                assert(clustnum == SummaryStruct.birds(1).neurons(nn).clustnum, 'asdfaf');
                
                
                MotifStats.neurons(nn).motif(j).SegmentsExtract = ...
                    lt_neural_SmoothFR(MotifStats.neurons(nn).motif(j).SegmentsExtract, ...
                    clustnum);
                
                
                %% ==== REMOVE ANY EXPERIMENTS THAT HAVE FR OF 0 DURING WINDOW OF INTEREST
                if isfield(MotifStats.neurons(nn).motif(j).SegmentsExtract, 'FRsmooth_xbin')
                    
                if RemoveTrialsZeroFR==1
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    tmp = segextract(1).FRsmooth_xbin;
                    premotorInds = find(tmp>(MotifStats.params.motif_predur + premotorWind(1)) ...
                        & tmp<(MotifStats.params.motif_predur + premotorWind(2)));
                    
                    alltrialFR = [segextract.FRsmooth_rate_CommonTrialDur];
                    alltrialFR = alltrialFR(premotorInds, :);
                    
                    trialstoremove = sum(alltrialFR,1)==0;
                    
                    if any(trialstoremove)
                        
%                         MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract(trialstoremove) = [];
                        MotifStats.neurons(nn).motif(j).SegmentsExtract(trialstoremove) = [];

                        disp(['REMOVED ' num2str(sum(trialstoremove)) '/' num2str(length(trialstoremove)) '(REASON: fr = 0)']);
                    end
                end
                end
                
                
                
            end
        end
        
        %%
        % === stick back into main structure
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS = MotifStats;
        
    end
end
disp('DONE! ---');


