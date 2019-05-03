function DirStruct = lt_DirSong_ExtractAFPbias(DirStruct)
%% lt 4/5/19 - collects average of day-means (each day subtract dir from u ndir)
% need to first run lt_DirSong_MotifID to extract DirStruct


%%

nbirds = length(DirStruct.bird);
for i=1:nbirds
   nmotif = length(DirStruct.bird(i).DAT.motifID);
   for mm=1:nmotif
      
       if isempty(DirStruct.bird(i).DAT.motifID(mm).rendnum)
           continue
       end
       ff = [DirStruct.bird(i).DAT.motifID(mm).rendnum.ff];
       tvals = [DirStruct.bird(i).DAT.motifID(mm).rendnum.datenum_song_SecRes];
       isdir = [DirStruct.bird(i).DAT.motifID(mm).rendnum.isDIR];
       
       ffmean_undir = grpstats(ff(isdir==0), floor(tvals(isdir==0)));
        ffmean_dir = grpstats(ff(isdir==1), floor(tvals(isdir==1)));
        
        afpbias = mean(ffmean_undir - ffmean_dir);
        
       DirStruct.bird(i).DAT.motifID(mm).means.afpbias = afpbias;
   end
end




