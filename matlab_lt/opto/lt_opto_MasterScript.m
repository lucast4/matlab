%% lt 10/11/18 - once do extractions for each bird, this plots summaries across birds

%% ============ LIST OF EXPERIEMNTS
ind =0;
clear ExptList;

% ================== 
ind=ind+1;
ExptList(ind).birdname = 'bu6bu98';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/bu6bu98/Opto_Stim_analy/Reversion1';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '29Sep2018-2021', 'of-up', ...
    '04Oct2018-2400', 'up-dn', ...
    '09Oct2018-2400', 'dn-up'...
    };

% ================== 
ind=ind+1;
ExptList(ind).birdname = 'pu83wh58';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/pu83wh58/Opto_Stim_analy/Reversion1';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '29Sep2018-2016', 'of-dn', ...
    '06Oct2018-2400', 'dn-up', ...
    '11Oct2018-2400', 'up-dn', ...
    };



%% ############################# PLOT INDIVIDUAL EXPERIEMNTS
close all;
expttoplot = 2;
OUTSTRUCT = lt_opto_ExtrBirdDat(ExptList(expttoplot).dirtoplot, ...
    ExptList(expttoplot).twind, ExptList(expttoplot).SwitchTimes);

%% ############################# SUMMARIZE ACROSS MULTIPLE EXPERIMENTS
% =================== PARAMS
normeachexpt=1; % NORMALIZE EACH EXPT (by mean of means over pos and
% neg train

% ===============
OUTSTRUCT = struct;

for j=1:length(ExptList)
    
    
 structbird = lt_opto_ExtrBirdDat(ExptList(j).dirtoplot, ...
    ExptList(j).twind, ExptList(j).SwitchTimes, ...
    0);
    
OUTSTRUCT.exptnum(j).dat=structbird;

end


% ================= PLOT, COMBINE ALL EXPT
allffstim = [];
alltraindir = [];
allexptnum = [];

for j=1:length(ExptList)
    
   ffdiff_stim = OUTSTRUCT.exptnum(j).dat.All_FFmeanStimMinusNostim;
   traindir = OUTSTRUCT.exptnum(j).dat.All_traindir;

   
% --------------------- NORMALIZE EACH EXPT (by mean of means over pos and
% neg train
if normeachexpt==1
    
    y = grpstats(ffdiff_stim(traindir~=0), traindir(traindir~=0), {'mean'});
    assert(length(y)==2, 'then doesnt have both upa nd dn train ...');
    normval = mean(y);
    ffdiff_stim = ffdiff_stim-normval;    
end

allffstim = [allffstim; ffdiff_stim];
alltraindir = [alltraindir; traindir];
allexptnum = [allexptnum; ones(size(traindir,1),1)*j];
      
end



lt_figure; hold on;
pcolexpt = lt_make_plot_colors(max(allexptnum),0,  0);
% ------
lt_subplot(2,2,1); hold on;
title('all dat');
scatter(alltraindir, allffstim, [], cell2mat(pcolexpt(allexptnum)'));
x = unique(alltraindir);
[y, ysem] = grpstats(allffstim, alltraindir, {'mean', 'sem'});
lt_plot(x+0.15, y, {'Errors', ysem, 'Color', 'k'});
xlim([-2 2]);
lt_plot_zeroline;


% ------
lt_subplot(2,2,2); hold on;
title('expt is dat');
ffall = [];
tdall =[];
xlabel('traindir');
ylabel('ff (stim - nostim)');
for j=unique(allexptnum')
    
    ffdiff_this = grpstats(allffstim(allexptnum==j), alltraindir(allexptnum==j), {'mean'});
    traindir_this = unique(alltraindir(allexptnum==j));
    
    plot(traindir_this, ffdiff_this, '-o', 'Color', pcolexpt{j});

    ffall = [ffall; ffdiff_this];
    tdall =[tdall; traindir_this];

end
[y, ysem] = grpstats(ffall, tdall, {'mean', 'sem'});
lt_plot(unique(tdall)+0.15, y, {'Errors', ysem});
lt_plot_zeroline;




