function MOTIFSTATS_Compiled = lt_neural_POPLEARN_ExtrWNtime(timewindhit, SummaryStruct,...
    MOTIFSTATS_Compiled, PARAMS)

%% extracts time of WN hit from raw data. LT 2018/2019


numbirds = length(MOTIFSTATS_Compiled.birds);

for i=1:numbirds
    
    numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    
    for ii=1:numneurons
        
        nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif);
        exptid = MOTIFSTATS_Compiled.SummaryStruct.birds(i).neurons(ii).exptID;
        
        hasdir = 0;
        for mm=1:nummotifs
            
            
            segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
            
            if isempty(segextract)
                continue
            end
            
            % ========== get time of WN for all trials
            ntrials = length(segextract);
            for tt = 1:ntrials
                fname = segextract(tt).song_filename;
                disp(fname);
                ons = segextract(tt).WithinSong_TokenOns;
                
                datdur = segextract(tt).actualmotifdur;
                
                if isempty(timewindhit)
                    windon = ons-PARAMS.motif_predur;
                    windoff = ons + datdur + PARAMS.motif_postdur;
                else
                    windon = ons + timewindhit(1);
                    windoff = ons + timewindhit(2);
                end
                [a] = fileparts(SummaryStruct.birds(i).neurons(ii).dirname);
                tmp = load([a '/' fname '.wntime.mat']);
                tmp.wnstruct;
                
                % ---- FIND WN ONSETS/OFFSETS WITHIN WINDOW (collect if on
                % OR off is within data window)
%                 assert(segextract(tt).global_ontime_motifInclFlank - (ons-PARAMS.motif_predur)<0.001, 'shouldnt this be the start of data window?');
                
                indsthis1 = tmp.wnstruct.WNonsets>windon & tmp.wnstruct.WNonsets<windoff;
%                 onsthis = tmp.wnstruct.WNonsets(indsthis1) - windon;
                onsthis = tmp.wnstruct.WNonsets(indsthis1) - (ons-PARAMS.motif_predur); % relative to start of the data window
                
                
                indsthis2 = tmp.wnstruct.WNoffsets>windon & tmp.wnstruct.WNoffsets<windoff;
%                 offthis = tmp.wnstruct.WNoffsets(indsthis2) - windon;
                offthis = tmp.wnstruct.WNoffsets(indsthis2) - (ons-PARAMS.motif_predur);
                
                wasTrialHit = any([indsthis1 indsthis2]);
                
                segextract(tt).hit_WN = wasTrialHit;
                segextract(tt).WNonset_sec = onsthis;
                segextract(tt).WNoffset_sec = offthis;
            end
            MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract = segextract;
        end
    end
end


