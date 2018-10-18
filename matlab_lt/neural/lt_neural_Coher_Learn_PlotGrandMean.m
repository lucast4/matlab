
%% lt 10/13/18 - each motif one val (mean over all chans), average over all expts
% ======= 1) EXTRACT
close all;
plotON=0;
averagechanpairs=1; % for each motif, average over all chan pairs
OUTSTRUCT = lt_neural_Coher_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS);

%% ======= 2) PLOT [ONE PLOT FOR EACH CASE (E.G. MOTIF/CHANPAIR])
close all;

useabs = 0; % if 1, then absolute values (wn minus diff)
plotON =0;
vssametype=0;
Yallswitch = lt_neural_Coher_Learn_GetSwitchCoh(SwitchStruct, OUTSTRUCT,...
    useabs, plotON, PARAMS, vssametype);


%% PLOT OVER ALL EXPERIMENTS
% TWO METHODS, 1) ALL CHAN PAIRS, OR 2) FOR EACH MOTIF FIRST AVERAGE OVER
% CHANNEL PAIRS ...
% ======================= PLOT MEANS
lt_figure; hold on;

% --------- 1) mean heat maps
cohmat = lt_neural_Coher_Cell2Mat(Yallswitch(:,1));

lt_subplot(3,2,1); hold on;
title('mean over switch');
ylabel('TARG(WN-base)');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);

lt_subplot(3,2,2); hold on;
% title('mean over switch');
ylabel('TARG');
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim);
lt_plot_zeroline;


% --------- 1) mean heat maps
cohmat = lt_neural_Coher_Cell2Mat(Yallswitch(:,2));

lt_subplot(3,2,3); hold on;
title('mean over switch');
ylabel('NONTARG(WN-base)');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);

lt_subplot(3,2,4); hold on;
% title('mean over switch');
ylabel('NONTARG');
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim);
lt_plot_zeroline;


% ------------
cohmat1 = lt_neural_Coher_Cell2Mat(Yallswitch(:,1)); 
cohmat2 = lt_neural_Coher_Cell2Mat(Yallswitch(:,2));
cohmat = cohmat1 - cohmat2;

lt_subplot(3,2,5); hold on;
title('mean over switch');
ylabel('TARG-NONTARG(WN-base)');
% lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);

lt_subplot(3,2,6); hold on;
% title('mean over switch');
ylabel('TARG-NONTARG(WN-base)');
% lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 2, '-', clim);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
lt_plot_zeroline;
% for j=1:size(cohmat,3)
%     cohmatthis = cohmat(:,:,j);
%         
%     
% end
