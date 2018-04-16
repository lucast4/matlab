function MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled)

%% lt 4/1/18 - removes all directed song trials from MOTIFSTATS_COMPILED

numbirds = length(MOTIFSTATS_Compiled.birds);
numdir = [];
numtot = [];
for i=1:numbirds
   
    numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
    
    for ii=1:numneurons
       
        nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif);
        
        for mm=1:nummotifs
           
            segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
            
            if isempty(segextract)
                continue
            end
            % ==== find inds of dir songs
            indsdir = [segextract.DirSong];
            disp(['DIR/TOT: ' num2str(sum(indsdir)) '/' num2str(length(indsdir))]);
            
            % ==== return to structure without dir
            MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract ...
                = segextract(~indsdir);
            
            numdir = [numdir sum(indsdir)];
            numtot = [numtot length(indsdir)];
        end
    end
end

% === add flag to indicate that DIR is removed
MOTIFSTATS_Compiled.DirSongRemoved = 'yes';
disp(['[ALL] DIR/TOT: ' num2str(sum(numdir)) '/' num2str(sum(numtot))]);

