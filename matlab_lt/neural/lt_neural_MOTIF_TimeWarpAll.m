function MOTIFSTATS_Compiled = lt_neural_MOTIF_TimeWarpAll(MOTIFSTATS_Compiled)

%% lt 4/11/18 - time warps from end to end all motifs in MOTIFSTATS_Compiled.
% assumes that extraction was from start to end of motif (i.e.
% Params_regexp.preAndPostDurRelSameTimept = 0);
% ALSO: time warps to match all neurons for a given motif.


if (0)
    %% ======================= LINEAR TIME WARP
    TimeWarpParams = {'pu69wh78', '(j)jjbhhg', [1:13], ...
        'pu69wh78', '(a)abhhg', [1:11]};
    NumBirds = length(MOTIFSTATS_Compiled.birds);
    
    MOTIFSTATS_Compiled.TimeWarpParams = TimeWarpParams;
    for i=1:NumBirds
        
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
        
        nneur = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
        for ii=1:nneur
            nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif);
            for iii=1:nummotifs
                
                motifthis = motiflist{iii};
                birdthis = MOTIFSTATS_Compiled.birds(i).birdname;
                
                ind1 = find(strcmp(TimeWarpParams, birdthis));
                ind2 = find(strcmp(TimeWarpParams, motifthis));
                
                ind3 = ind1(ind1 == ind2-1)+2; % actual ind of params
                
                if isempty(ind3)
                    % then this bird or motif is not specificed in params
                    disp(['PROBLEM - b ' birdname '- motif ' num2str(motifthis) ' NOT SPECIFIED (WILL NOT TIME WARP)']);
                    continue
                end
                
                disp([birdthis '-n' num2str(ii) '-' motifthis]);
                
                % ================= DO TIME WARP
                regionstowarp = TimeWarpParams{ind3};
                segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract;
                
                expectedsegs = 2*length(segextract(1).matchlabel) - 1;
                segextract = lt_neural_LinTimeWarpSegmented(segextract, ...
                    regionstowarp, expectedsegs);
                
                MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract = ...
                    segextract;
            end
        end
    end
    
else
    %% ==================== LINEAR TIME WARP, COMBINGIN ALL NEURONS
    % ================== FIRST, CHECK PARAMS
    assert(MOTIFSTATS_Compiled.birds(1).Params_regexp.preAndPostDurRelSameTimept==0, 'assumes have entire motif data...');
    
    
    % =========================== RUN
%     TimeWarpParams = {'pu69wh78', '(j)jjbhhg', [1:13], ...
%         'pu69wh78', '(a)abhhg', [1:11], ...
%         'wh44wh39', '(n)hh', [1:5], ...
%         'wh44wh39', '(d)kccbb', [1:11]};
    NumBirds = length(MOTIFSTATS_Compiled.birds);
    
%     MOTIFSTATS_Compiled.TimeWarpParams = TimeWarpParams;
    
    for i=1:NumBirds
        
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
        nummotifs = length(motiflist);
        nneur = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
        
        for iii=1:nummotifs
            
            motifthis = motiflist{iii};
            birdthis = MOTIFSTATS_Compiled.birds(i).birdname;

            % ========== figure out what regions to warp based on motif
%             indpar = strfind(motifthis, ')');
%             nsyls = length(motifthis) - indpar +1
            nsyls = length(motifthis)-2; % -2 to take into account parantheses
            regionstowarp = 1:(nsyls + nsyls-1); % i.e all syls and gaps
            
%             ind1 = find(strcmp(TimeWarpParams, birdthis));
%             ind2 = find(strcmp(TimeWarpParams, motifthis));
%             ind3 = ind1(ind1 == ind2-1)+2; % actual ind of params
%             
%             if isempty(ind3)
%                 % then this bird or motif is not specificed in params
%                 disp(['PROBLEM - b ' birdname '- motif ' num2str(motifthis) ' NOT SPECIFIED (WILL NOT TIME WARP)']);
%                 continue
%             end

            % ================== combine segextract across all neurons
            segextractAll = [];
            segIndsAll = cell(1, nneur);  % --- save inds to be able to put back into specific neurons
            matchlabel = '';
            for ii=1:nneur
                disp([birdthis '-n' num2str(ii) '-' motifthis]);
                
                % =================
                segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract;
                if isempty(segextract)
                    continue
                end
                segIndsAll{ii} = length(segextractAll)+1:(length(segextractAll)+length(segextract));
                
                % --- add to concated
                segextractAll = [segextractAll segextract];
                
                matchlabel = segextract(1).matchlabel;
            end
            
            % ============== DO TIMEWARP
%             regionstowarp = TimeWarpParams{ind3};
            
            expectedsegs = 2*length(matchlabel) - 1;

            segextractAll = lt_neural_LinTimeWarpSegmented(segextractAll, ...
                regionstowarp, expectedsegs);
            
            % =========== SLIDE BACK INTO MOTIFSTATS
            for ii=1:nneur
                
                indstotake = segIndsAll{ii};
                assert(length(indstotake) == length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract), 'asfsda');
                
                MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(iii).SegmentsExtract = ...
                    segextractAll(indstotake);
            end
            
            
        end
    end
end