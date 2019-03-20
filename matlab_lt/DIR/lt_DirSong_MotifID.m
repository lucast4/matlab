function dirstruct = lt_DirSong_MotifID
%% lt 3/14/19 - gets unique motifID,trial by trial data (written fgor nerual experiments)

dirstruct = lt_DirSong_Extract;

%%
nbirds = length(dirstruct.bird);
gethitsyls = 0;
for i=1:nbirds
   bname = dirstruct.bird(i).birdname;
   
   [~, ~, ~, motiflist] = lt_neural_QUICK_MotifID(bname);
   
   % ======= for all motifs, extract ff
   ListOfDirs_UNDIR = dirstruct.bird(i).UNDIR_directories;
   ListOfDirs_DIR= dirstruct.bird(i).DIR_directories;
   ListOfBatch = cell(1, 2*length(ListOfDirs_UNDIR));
   ListOfBatch(:)= {'batchall.LABELED'};
   MotifsToExtract = {};
   for j=1:length(motiflist)
%        motiflist{j};
       MotifsToExtract = [MotifsToExtract motiflist{j}];
   end
%    MotifsToExtract = motiflist;
   
    DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_UNDIR, ListOfDirs_DIR, ...
        ListOfBatch, MotifsToExtract, gethitsyls);
    
    % =============== ASSIGN TRIAL VALUES FOR EACH UNIQUE MOTIF
    for j=1:length(motiflist)
       
%         indsthis = strcmp({DATSTRUCT.motif.motif}, motiflist{j});
        indsthis = find(ismember({DATSTRUCT.motif.motif}, motiflist{j}));
        
        tmp = [DATSTRUCT.motif(indsthis).rendnum]; % concatenate all mo0tifs.
        if isempty(tmp)
            continue
        end
        
        % --- get only unique trials, some motifs probably overlap
        [indsgrp] =  lt_tools_grp2idx({[tmp.datenum_song_SecRes]', ...
            [tmp.time_withinsong]'}); % each trial gets an ID.
        [~, indsgood] = unique(indsgrp); % get indices for the unqiue trials

        datthis = tmp(indsgood);
        
        % ==================== SAVE OUTPUT
        dirstruct.bird(i).DAT.motifID(j).rendnum = datthis;
    end
end


