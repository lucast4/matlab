% plot distribution of max repeats during song bouts
% function [maxReps] = ek_syllDist(syll, batchf)

syll = 'a';
syll2 = 'b';
batchf = 'batch.labeled.all';

fid = fopen(batchf);
fline = fgetl(fid);

maxReps = [];
reps = [];
while ischar(fline)

   disp(fline);
   
   datnotmat = load([fline '.not.mat']);
   
%    inds_a = strfind(datnotmat.labels, strtoget_a);
   song = strrep(datnotmat.labels, 'j', '-');
   C = strsplit(song, '-');
   nReps = [];
   for i = 1 : length(C)
       B = strsplit(C{i}, syll2);
       for j = 1 : length(B)
           idx = strfind(B{j}, syll);
           if isempty(idx) == 0
               nReps = [nReps length(idx)];
           end
       end
   end
    maxReps = [maxReps max(nReps)];
    reps = [reps nReps];
    fline = fgetl(fid);
end

% plot histogram
figure(1)
hold on
% nhist(maxReps, 'color', 'white')
nhist(maxReps, 'color', 'summer')
title(strcat('Max Number of Repeats per Bout: ', syll))
xlabel('Max Repeats')
ylabel('Number of Song Bouts')
set(gca,'fontsize', 12, 'TickDir','out');
% 
% % maxVal = mode(reps);
% % peak = length(find(reps == maxVal));
% figure(2)
% hold on
% % nhist(reps, 'color', 'white')
% nhist(reps, 'color', 'summer')
% title(strcat('Repeats per Motif: ', syll))
% xlabel('# Reps')
% ylabel('Number of Motifs')
% set(gca,'fontsize', 12, 'TickDir','out');

% %%
% figure(2) % normalized probability distribution of repeats per motif
% hold on
% % nhist(reps, 'color', 'white')
% % nhist(repDist, 'color', 'qualitative')
% nhist(repDist,'samebins', 'color', 'qualitative', 'proportion', 1) 
% % lt_plot_histogram(reps)
% title(strcat('Repeats per Motif: ', syll))
% xlabel('# Reps')
% ylabel('Probability')
% set(gca,'fontsize', 12, 'TickDir','out');
% 
% figure(3) % normalized probability distribution of repeats per motif
% hold on
% % nhist(reps, 'color', 'white')
% % nhist(repDist, 'color', 'qualitative')
% nhist(maxRepDist,'samebins', 'color', 'qualitative', 'proportion', 1) 
% % lt_plot_histogram(reps)
% title(strcat('Max Repeats per Motif: ', syll))
% xlabel('# Reps')
% ylabel('Probability')
% set(gca,'fontsize', 12, 'TickDir','out');











