%% lt 10/11/18 - once do extractions for each bird, this plots summaries across birds

%% ============ LIST OF EXPERIEMNTS
clear all; close all;
ind = 0;

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
ExptList(ind).TrainDirForStimContext = {...
    };


% ================== 
ind=ind+1;
ExptList(ind).birdname = 'bu6bu98';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/bu6bu98/Opto_Stim_analy/Reversion2';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '14Nov2018-1410', 'of-dn', ...
    '19Nov2018-1339', 'dn-up', ...
    };
ExptList(ind).TrainDirForStimContext = {...
    };


% ================== 
ind=ind+1;
ExptList(ind).birdname = 'pu83wh58v1';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/pu83wh58/Opto_Stim_analy/Reversion1';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '29Sep2018-2016', 'of-dn', ...
    '06Oct2018-2400', 'dn-up', ...
    '11Oct2018-2400', 'up-dn', ...
    };
ExptList(ind).TrainDirForStimContext = {...
    };


% ================== 
ind=ind+1;
ExptList(ind).birdname = 'pu83wh58v2';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/pu83wh58/Opto_Stim_analy/OneDirLearn';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '07Nov2018-2020', 'of-up', ...
    '14Nov2018-2330', 'up-dn', ...
    '20Nov2018-1405', 'dn-up', ...
    };
ExptList(ind).TrainDirForStimContext = {...
    };


% ######################################### [OLD EXPT]
ind=ind+1;
ExptList(ind).birdname = 'or60';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/or60/Opto_Stim_analy/Reversion1';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '19Feb2017-1115', 'of-dn', ...
    '22Feb2017-1737', 'dn-up',...
    '06Mar2017-1137', 'up-dn',...
    '09Mar2017-1223', 'dn-up', ...
    '09Mar2017-1821', 'up-dn',...
    '10Mar2017-1023', 'dn-up',...
    '10Mar2017-1338', 'up-of'
    };
ExptList(ind).TrainDirForStimContext = {...
    };


ind=ind+1;
ExptList(ind).birdname = 'or60';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/or60/Opto_Stim_analy/Reversion2';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '13Mar2017-2312', 'of-up', ...
    '15Mar2017-2300', 'up-dn', ...
    '18Mar2017-1432', 'dn-up'};
ExptList(ind).TrainDirForStimContext = {...
    };


% =========== TO ADD: WH73, 
ind=ind+1;
ExptList(ind).birdname = 'wh73pk60';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/wh73pk61/Opto_Stim_analy/Reversion';
ExptList(ind).twind = 4;
ExptList(ind).SwitchTimes = {...
    '08Mar2015-2330', 'of-up', ...
    '25Mar2015-1832', 'up-of', ... % IMPORTANT: from before and including 25Mar2015-1832 I did not enter all swicfhes. many switches, each day, and multiple syls sometimes in diff directions. fine as I am not looking at that data
    '26Mar2015-1233', 'of-up', ...
    '28Mar2015-1714', 'up-dn', ...
    '28Mar2015-1920', 'dn-up', ...
    '29Mar2015-1551', 'up-dn', ...
    '29Mar2015-1838', 'dn-up', ...
    '30Mar2015-1650', 'up-dn', ...
    '30Mar2015-1856', 'dn-up', ...
    '31Mar2015-1122', 'up-dn', ...
    '31Mar2015-1550', 'dn-up', ...
    '01Apr2015-1238', 'up-dn', ...
    '01Apr2015-1548', 'dn-up', ...
    '01Apr2015-1925', 'up-dn', ...
    '01Apr2015-2330', 'dn-up', ...
    '02Apr2015-1143', 'up-of', ... % NOTE: actually started assoc here... [and other stuff...]
    '20Apr2015-1232', 'of-dn', ...
    '01May2015-1434', 'dn-of', ...
    '06May2015-1131', 'of-up', ...
    '14May2015-2330', 'up-dn', ...
    '16May2015-1946', 'dn-of', ...
   };
ExptList(ind).TrainDirForStimContext = {...
    };


%% ============== LIST OF EXPERIMENTS [ASSOCIATION]
% NOTE: for these experiments make switchtime mean the direction of
% learning that laser is associated with.
clear all; close all;
ind =0;

% ==================
ind=ind+1;
ExptList(ind).birdname = 'bu6bu98';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/bu6bu98/Opto_Stim_analy/Association2';
ExptList(ind).twind = 1;
% ExptList(ind).SwitchTimes = {...
%     '12Oct2018-1833', 'of-dn' ...
%     '20Oct2018-2313', 'dn-up' ...
%     '27Oct2018-1535', 'up-dn' ...
%     '07Nov2018-0914', 'dn-up' ...
%     '13Nov2018-2330', 'up-dn' ...    
%     '14Nov2018-1410', 'dn-of' ...    % ACTUALLY TURNED OFF LASER...
% };
% ExptList(ind).daystoignore = {...
%     '04Nov2018', ...
%     '05Nov2018', ...
%     '06Nov2018', ...
%     '07Nov2018', ...
%     }; % DAYS WHEN STIM IS OFF
ExptList(ind).SwitchTimes = {...
    '12Oct2018-1833', 'of-dn' ...
    '20Oct2018-2313', 'dn-up' ...
    '27Oct2018-1535', 'up-dn' ...
    '03Nov2018-2247', 'dn-of' ...
    '07Nov2018-2014', 'of-up' ...
    '13Nov2018-2330', 'up-dn' ...    
    '14Nov2018-1410', 'dn-of' ...    % ACTUALLY TURNED OFF LASER...
}; % NOTE: this now indicates both STIM changes (ON/OFF) and WN direction changes...
ExptList(ind).daystoignore = {...
    }; % DAYS WHEN STIM IS OFF



% ==================
ind=ind+1;
ExptList(ind).birdname = 'pu83wh58v1';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/pu83wh58/Opto_Stim_analy/Association2';
ExptList(ind).twind = 1;
% ExptList(ind).SwitchTimes = {...
%     '15Oct2018-2330', 'of-up' ...
%     '21Oct2018-2036', 'up-dn' ...
%     '27Oct2018-1539', 'dn-up' ...
% };
% ExptList(ind).daystoignore = {...
%     '04Nov2018', ...
%     '05Nov2018', ...
%     '06Nov2018', ...
%     '07Nov2018', ...
%     }; % DAYS WHEN STIM IS OFF
ExptList(ind).SwitchTimes = {...
    '15Oct2018-2330', 'of-up' ...
    '21Oct2018-2036', 'up-dn' ...
    '27Oct2018-1539', 'dn-up' ...
    '03Nov2018-2247', 'up-of' ...
};
ExptList(ind).daystoignore = {...
    }; % DAYS WHEN STIM IS OFF





% ====================== [IGNORE, OLD VERSION BEFORE CLEANED UP LABELING]
% % PREVIOUSLY: had issues in that noise syllables in labeling. Cleaned up by
% % going thru and checking all autolabels (see below)
% ind=ind+1;
% ExptList(ind).birdname = 'wh73pk61';
% ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/wh73pk61/Opto_Stim_analy/Association';
% ExptList(ind).twind = 4;
% % ExptList(ind).SwitchTimes = {...
% %     '02Apr2015-1143', 'of-up', ...
% %     '04Apr2015-1017', 'up-dn', ...
% %     '06Apr2015-1241', 'dn-up', ...
% %     '10Apr2015-1837', 'up-of', ...
% %     '12Apr2015-2330', 'of-dn', ...
% %     '14Apr2015-2330', 'dn-up', ...
% %     '15Apr2015-2330', 'up-of', ...
% %     };
% % ExptList(ind).daystoignore = {...
% %     '11Apr2015', ...
% %     '12Apr2015', ...
% %     }; % DAYS WHEN STIM IS OFF
% ExptList(ind).SwitchTimes = {...
%     '02Apr2015-1143', 'of-up', ...
%     '04Apr2015-1017', 'up-dn', ...
%     '06Apr2015-1241', 'dn-up', ...
%     '10Apr2015-1837', 'up-of', ...
%     '12Apr2015-2330', 'of-dn', ...
%     '14Apr2015-2330', 'dn-up', ...
%     '15Apr2015-2330', 'up-of', ...
%     };
% ExptList(ind).daystoignore = {...
%     }; % DAYS WHEN STIM IS OFF



% ====================== [CLEANED UP - GOOD]
% [IGNORE, FOR FEW REASONS. 1) start with opnly WN on 50%. 2) then few days
% with almost not WN (did not update WN) (4/4 and 4/5); 3) many days with
% also "probe" epochs, which reduce the assoc with laser and WN. Even if
% include, stil only one(?) day with down(up?) association. Therefore is
% incredibly noisy and not enough data.
ind=ind+1;
ExptList(ind).birdname = 'wh73pk61';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/wh73pk61/Opto_Stim_analy/Association1';
ExptList(ind).twind = 1;
% ExptList(ind).SwitchTimes = {...
%     '02Apr2015-1143', 'of-up', ...
%     '04Apr2015-1017', 'up-dn', ...
%     '06Apr2015-1241', 'dn-up', ...
%     '10Apr2015-1837', 'up-of', ...
%     '12Apr2015-2330', 'of-dn', ...
%     '14Apr2015-2330', 'dn-up', ...
%     '15Apr2015-2330', 'up-of', ...
%     };
% ExptList(ind).daystoignore = {...
%     '11Apr2015', ...
%     '12Apr2015', ...
%     }; % DAYS WHEN STIM IS OFF
ExptList(ind).SwitchTimes = {...
    '02Apr2015-1143', 'of-up', ...
    '04Apr2015-1017', 'up-dn', ...
    '06Apr2015-1241', 'dn-up', ...
    '10Apr2015-1837', 'up-of', ...
    '12Apr2015-2330', 'of-dn', ...
    '14Apr2015-2330', 'dn-up', ...
    '15Apr2015-2330', 'up-of', ...
    };
ExptList(ind).daystoignore = {...
    }; % DAYS WHEN STIM IS OFF


% #########################################################
% ################### LOAD DATA EXTRACTED USING CONTEXT ANALYSIS
% SEE or60_analysis... for example.
% there, saves summary, ff difference and train dir for each day. Here
% extracts that.

ind=ind+1;
ExptList(ind).birdname = 'or60';
ExptList(ind).exptID = 'Association2';
ExptList(ind).ctxtsummarydat = '/bluejay5/lucas/analyses/opto/ConvertFromContextToOpto/dat_or60_Association2.mat';


%% ################## SETS OF EXPERIMENTS THAT HAVE BOTH REVERSION AND ASSOC.
% In alternating order (reversion, assoc, reverison ...)

clear all; close all;
ind = 0;
isrev = [1 1 1 1 1 0 0 0]'; % 1s and 0s, 1 for reverseion 0 for not

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
ExptList(ind).TrainDirForStimContext = {...
    };


% ================== 
ind=ind+1;
ExptList(ind).birdname = 'bu6bu98';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/bu6bu98/Opto_Stim_analy/Reversion2';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '14Nov2018-1410', 'of-dn', ...
    '19Nov2018-1339', 'dn-up', ...
    };
ExptList(ind).TrainDirForStimContext = {...
    };


% ================== 
ind=ind+1;
ExptList(ind).birdname = 'pu83wh58v1';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/pu83wh58/Opto_Stim_analy/Reversion1';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '29Sep2018-2016', 'of-dn', ...
    '06Oct2018-2400', 'dn-up', ...
    '11Oct2018-2400', 'up-dn', ...
    };
ExptList(ind).TrainDirForStimContext = {...
    };

ind=ind+1;
ExptList(ind).birdname = 'or60';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/or60/Opto_Stim_analy/Reversion1';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '19Feb2017-1115', 'of-dn', ...
    '22Feb2017-1737', 'dn-up',...
    '06Mar2017-1137', 'up-dn',...
    '09Mar2017-1223', 'dn-up', ...
    '09Mar2017-1821', 'up-dn',...
    '10Mar2017-1023', 'dn-up',...
    '10Mar2017-1338', 'up-of'
    };
ExptList(ind).TrainDirForStimContext = {...
    };


ind=ind+1;
ExptList(ind).birdname = 'or60';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/or60/Opto_Stim_analy/Reversion2';
ExptList(ind).twind = 1;
ExptList(ind).SwitchTimes = {...
    '13Mar2017-2312', 'of-up', ...
    '15Mar2017-2300', 'up-dn', ...
    '18Mar2017-1432', 'dn-up'};
ExptList(ind).TrainDirForStimContext = {...
    };


% ==================
ind=ind+1;
ExptList(ind).birdname = 'bu6bu98';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/bu6bu98/Opto_Stim_analy/Association2';
ExptList(ind).twind = 1;
% ExptList(ind).SwitchTimes = {...
%     '12Oct2018-1833', 'of-dn' ...
%     '20Oct2018-2313', 'dn-up' ...
%     '27Oct2018-1535', 'up-dn' ...
%     '07Nov2018-0914', 'dn-up' ...
%     '13Nov2018-2330', 'up-dn' ...    
%     '14Nov2018-1410', 'dn-of' ...    % ACTUALLY TURNED OFF LASER...
% };
% ExptList(ind).daystoignore = {...
%     '04Nov2018', ...
%     '05Nov2018', ...
%     '06Nov2018', ...
%     '07Nov2018', ...
%     }; % DAYS WHEN STIM IS OFF
ExptList(ind).SwitchTimes = {...
    '12Oct2018-1833', 'of-dn' ...
    '20Oct2018-2313', 'dn-up' ...
    '27Oct2018-1535', 'up-dn' ...
    '03Nov2018-2247', 'dn-of' ...
    '07Nov2018-2014', 'of-up' ...
    '13Nov2018-2330', 'up-dn' ...    
    '14Nov2018-1410', 'dn-of' ...    % ACTUALLY TURNED OFF LASER...
}; % NOTE: this now indicates both STIM changes (ON/OFF) and WN direction changes...
ExptList(ind).daystoignore = {...
    }; % DAYS WHEN STIM IS OFF



% ==================
ind=ind+1;
ExptList(ind).birdname = 'pu83wh58v1';
ExptList(ind).dirtoplot = '/bluejay5/lucas/birds/pu83wh58/Opto_Stim_analy/Association2';
ExptList(ind).twind = 1;
% ExptList(ind).SwitchTimes = {...
%     '15Oct2018-2330', 'of-up' ...
%     '21Oct2018-2036', 'up-dn' ...
%     '27Oct2018-1539', 'dn-up' ...
% };
% ExptList(ind).daystoignore = {...
%     '04Nov2018', ...
%     '05Nov2018', ...
%     '06Nov2018', ...
%     '07Nov2018', ...
%     }; % DAYS WHEN STIM IS OFF
ExptList(ind).SwitchTimes = {...
    '15Oct2018-2330', 'of-up' ...
    '21Oct2018-2036', 'up-dn' ...
    '27Oct2018-1539', 'dn-up' ...
    '03Nov2018-2247', 'up-of' ...
};
ExptList(ind).daystoignore = {...
    }; % DAYS WHEN STIM IS OFF


% #########################################################
% ################### LOAD DATA EXTRACTED USING CONTEXT ANALYSIS
% SEE or60_analysis... for example.
% there, saves summary, ff difference and train dir for each day. Here
% extracts that.

ind=ind+1;
ExptList(ind).birdname = 'or60';
ExptList(ind).exptID = 'Association2';
ExptList(ind).ctxtsummarydat = '/bluejay5/lucas/analyses/opto/ConvertFromContextToOpto/dat_or60_Association2.mat';

%% ################# MAKE UNIQUE ID FOR EACH EXPT
for i=1:length(ExptList)
    if ~isempty(ExptList(i).dirtoplot)
   [~, tag] = fileparts(ExptList(i).dirtoplot);
    end
    ExptList(i).exptID = tag;
end

%% ======================= [DISPLAY] Show list of all birds/experiments
disp('-------------------');
for j=1:length(ExptList)
   disp([num2str(j) ' == ' ExptList(j).birdname  '-' ExptList(j).exptID]);
end

%% ############################# PLOT INDIVIDUAL EXPERIEMNTS

close all;
expttoplot = 2;
onlylongepoch = 0;
StartDaySkipTime = 1;
lt_opto_ExtrBirdDat(ExptList(expttoplot).dirtoplot, ...
    ExptList(expttoplot).twind, ExptList(expttoplot).SwitchTimes, ...
    1, StartDaySkipTime, onlylongepoch);
title([ExptList(expttoplot).birdname '-' ExptList(expttoplot).exptID]);


%% ############################# SUMMARIZE ACROSS MULTIPLE EXPERIMENTS
%% ============= [EXTRACTION]
% =================== PARAMS
normeachexpt=0; % NORMALIZE EACH EXPT (by mean of means over pos and
% neg train
valfield = 'All_FFmeanStimMinusNostim';
% valfield = 'All_Dprime';
% valfield = 'All_FFmedianStimMinusNostim';
groupbybird = 0; % if 0, then by expt.
onlylongepoch = 0;
StartDaySkipTime = 1;

[allffstim, alltraindir, allexptnum, allbirdname] ...
    = lt_opto_ExtrAllDat(ExptList, normeachexpt, valfield, groupbybird, ...
    onlylongepoch, StartDaySkipTime);


%% ================ [PLOTS] Summary, each bird histrogram of stim effects.
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

numexpts = max(allexptnum);
xcenters = min(allffstim)-1:2:max(allffstim)+1;
for j=1:numexpts
   
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

if groupbybird==1
    bname = unique(allbirdname(allexptnum==j));
    plottit = bname{1};
else
    plottit = [ExptList(j).birdname '-' ExptList(j).exptID];
end
title(plottit);
xlabel('ff(stim-nostim) [bu=UP days]');


ff = allffstim(allexptnum==j);
td = alltraindir(allexptnum==j);
% -- train up
pcol = 'b';
lt_plot_histogram(ff(td==1), xcenters, 1, 0, 0.5, 1, pcol);

% -- train dn
pcol = 'r';
lt_plot_histogram(ff(td==-1), xcenters, 1, 0, 0.5, 1, pcol);
    

end



lt_figure; hold on;

numexpts = max(allexptnum);
xcenters = min(allffstim)-1:2:max(allffstim)+1;

Y = {};
X = [];
for j=1:numexpts
    
    
    if groupbybird==1
        bname = unique(allbirdname(allexptnum==j));
        plottit = bname{1};
    else
        plottit = [ExptList(j).birdname '-' ExptList(j).exptID];
    end
    lt_plot_text(j*2-0.5, 35, plottit);
    xlabel('ff(stim-nostim) [bu=UP days]');
    
    ff = allffstim(allexptnum==j);
    td = alltraindir(allexptnum==j);
    
    % -- train dn
    Y = [Y; ff(td==-1)];
    
    % -- train up
    Y = [Y; ff(td==1)];
    
    X = [X; [j*2-1+0.3 j*2-0.3]'];
end


Ymean = cellfun(@mean, Y);
Ysem = cellfun(@lt_sem, Y);
% -- dn
lt_plot_bar(X(1:2:end), Ymean(1:2:end), {'Errors', Ysem(1:2:end), 'Color', 'r', ...
    'BarWidth', 0.25});
% -- up
lt_plot_bar(X(2:2:end), Ymean(2:2:end), {'Errors', Ysem(2:2:end), 'Color', 'b', ...
    'BarWidth', 0.25});
lt_plot_MultDist(Y, X, 0, 'k', 1);
lt_plot_zeroline;
%% ================ [PLOTS] Summary plots across experiments
lt_figure; hold on;
pcolexpt = lt_make_plot_colors(max(allexptnum),0,  0);
% ------
lt_subplot(2,2,1); hold on;
title('all dat');
try
scatter(alltraindir, allffstim, [], cell2mat(pcolexpt(allexptnum)'));
catch err
scatter(alltraindir, allffstim, [], cell2mat(pcolexpt(allexptnum)));
end
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


% ------
lt_subplot(2,2,3); hold on;
title('expt is dat');
ffall = [];
tdall =[];
xlabel('traindir');
ylabel('ff (stim - nostim)');
for j=unique(allexptnum')
    
    ffdiff_this = grpstats(allffstim(allexptnum==j), alltraindir(allexptnum==j), {'mean'});
    traindir_this = unique(alltraindir(allexptnum==j));
    
    indtmp = traindir_this==0;
    ffdiff_this(indtmp) = [];
    traindir_this(indtmp) = [];
    
    plot(traindir_this, ffdiff_this, '-o', 'Color', pcolexpt{j});

    % --- plot name
    if groupbybird==0
        lt_plot_text(1.2, ffdiff_this(traindir_this==1), [ExptList(j).birdname '-' ExptList(j).exptID], pcolexpt{j});
    else
        bnamethis = unique(allbirdname(allexptnum==j));
        lt_plot_text(1.2, ffdiff_this(traindir_this==1), [bnamethis{1}], pcolexpt{j});
    end
    ffall = [ffall; ffdiff_this];
    tdall =[tdall; traindir_this];

end
[y, ysem] = grpstats(ffall, tdall, {'mean', 'sem'});
lt_plot(unique(tdall)+0.15, y, {'Errors', ysem});
lt_plot_zeroline;
[~, p] = ttest(ffall(tdall==-1), ffall(tdall==1));
lt_plot_pvalue(p, 'ttest, means', 1);


%% ============== [COMPARING REVERSION AND ASSOCIATION WITHIN SAME EXPERIMENT]

% ====================== EXTRACT [REVERSION]
normeachexpt=1; % NORMALIZE EACH EXPT (by mean of means over pos and
% neg train
valfield = 'All_FFmeanStimMinusNostim';
% valfield = 'All_Dprime';
% valfield = 'All_FFmedianStimMinusNostim';
groupbybird = 1; % if 0, then by expt.
onlylongepoch = 1;
StartDaySkipTime = 1;

[allffstim_REV, alltraindir_REV, allexptnum_REV, allbirdname_REV] ...
    = lt_opto_ExtrAllDat(ExptList(isrev==1), normeachexpt, valfield, groupbybird, ...
    onlylongepoch, StartDaySkipTime);


% ====================== EXTRACT [ASSOC]
normeachexpt=1; % NORMALIZE EACH EXPT (by mean of means over pos and
% neg train
valfield = 'All_FFmeanStimMinusNostim';
% valfield = 'All_Dprime';
% valfield = 'All_FFmedianStimMinusNostim';
groupbybird = 1; % if 0, then by expt.
onlylongepoch = 0;
StartDaySkipTime = 1;

[allffstim_ASSOC, alltraindir_ASSOC, allexptnum_ASSOC, allbirdname_ASSOC] ...
    = lt_opto_ExtrAllDat(ExptList(isrev==0), normeachexpt, valfield, groupbybird, ...
    onlylongepoch, StartDaySkipTime);

%% ======= [PLOT]

usedprime = 0;

% ============== SCATTER PLOT, EACH EXPERIMENT
lt_figure; hold on;
xlabel('reversion effect of stim [UPtrain - DNtrain]');
ylabel('association effect of stim [pairedUP - pairedDN]');
listofbirds = unique([allbirdname_REV; allbirdname_ASSOC]);
if usedprime==1
    title('dprime');
else
    title('diff of mean ff');
end


for j=1:length(listofbirds)
   birdthis = listofbirds{j};
   
   Y = [];
   Yerr = [];
   
   % ===== REVERSION
   indthis = strcmp(allbirdname_REV, birdthis);
   ff = allffstim_REV(indthis);
   tdir = alltraindir_REV(indthis);
  
   if usedprime==0
          Y = [Y; mean(ff(tdir==1)) - mean(ff(tdir==-1))];
   else
   Y = [Y; lt_tools_dprime(ff(tdir==1), ff(tdir==-1))];
   end
   
   % ===== REVERSION
   indthis = strcmp(allbirdname_ASSOC, birdthis);
   ff = allffstim_ASSOC(indthis);
   tdir = alltraindir_ASSOC(indthis);
  
   if usedprime==0
          Y = [Y; mean(ff(tdir==1)) - mean(ff(tdir==-1))];
   else
   Y = [Y; lt_tools_dprime(ff(tdir==1), ff(tdir==-1))];
   end

   % ======== plot
   lt_plot(Y(1), Y(2), {'Color', 'r'});
   
   % --- plot text of expt 
   lt_plot_text(Y(1), Y(2), birdthis, 'm');
end

lt_plot_zeroline;
lt_plot_zeroline_vert;


%% ################### [DOES EFFECT OF STIM CORRELATE WITH ACROSS DAY CHANGE?]
% ============ DO THIS WITH ASSOCIATION DATA
% NOTE:
% stim effects collect one per day. for or60 (context analysis) normalyl
% gets multipel balues per day, (isnce based on phases) so modifed here to
% fix that.
% confirmed that day means and stim effects data are all aligned (ie.. day
% by day correposntds).


% ============ 1) GET EFFECT OF STIM
normeachexpt=0; % NORMALIZE EACH EXPT (by mean of means over pos and
% neg train
valfield = 'All_FFmeanStimMinusNostim';
% valfield = 'All_Dprime';
% valfield = 'All_FFmedianStimMinusNostim';
groupbybird = 0; % if 0, then by expt.
onlylongepoch = 0;
StartDaySkipTime = 1;

[allffstim, ~, allexptnum, allbirdname, alldaynums] ...
    = lt_opto_ExtrAllDat(ExptList, normeachexpt, valfield, groupbybird, ...
    onlylongepoch, StartDaySkipTime);

% ==== get one value per day
indsgrp = lt_tools_grp2idx({allexptnum, allbirdname, alldaynums});
allffstim = grpstats(allffstim, indsgrp, 'mean');
allexptnum = grpstats(allexptnum, indsgrp, 'mean');
alldaynums = grpstats(alldaynums, indsgrp, 'mean');



% ============ 2) GET MEAN FF DAY BY DAY
OUTSTRUCT = struct;

SlopesAll_ffvsday = []; % [mean, CIl CIu]
StimEffectMean = []; % [mean, CIl CIu]

DayStats_StimEffect_MeanFF_DayNum = {}; % each cell one expt. each expt gest [StimEffect MeanFF DayNum], which is size ndays x 3

for j=1:length(ExptList)
    
    if isempty(ExptList(j).dirtoplot)
        % then tyr to extract previously saved summary data
        
        tmp = load(ExptList(j).ctxtsummarydat);
        structbird = tmp.dat;
        
        % ================== get day by day change in pitch
        ffall = [];
        tvalsall = [];
        for nn = 1:length(structbird.SORTED_DATA.ByNoteGroup)
            goodphases = find(~isnan(tmp.dat.PARAMS_GLOB.Assoc_PhaseLearndirMapping));

            %            ffvals = structbird.SORTED_DATA.ByNoteGroup(nn).Stats_OneDataPtPerEpoch.RAW.FFvals
            ffvals = cell2mat(cellfun(@transpose, structbird.SORTED_DATA.ByNoteGroup(nn).Stats_OneDataPtPerEpoch.RAW.FFvals, 'UniformOutput', 0));
            tvals = cell2mat(cellfun(@transpose, structbird.SORTED_DATA.ByNoteGroup(nn).Stats_OneDataPtPerEpoch.RAW.Tvals, 'UniformOutput', 0));
            phases = cell2mat(cellfun(@transpose, structbird.SORTED_DATA.ByNoteGroup(nn).Stats_OneDataPtPerEpoch.RAW.PhaseNumvals, 'UniformOutput', 0));
            
            % --- only keep those within good phaes - i.e phases inw hich
            % both WN and laser were on
            ffall = [ffall ffvals(ismember(phases, goodphases))];
            tvalsall = [tvalsall tvals(ismember(phases, goodphases))];
        end
        ffmeans_byday = grpstats(ffall, floor(tvalsall), {'mean'});
        days = unique(floor(tvalsall))';
    else
        structbird = lt_opto_ExtrBirdDat(ExptList(j).dirtoplot, ...
            ExptList(j).twind, ExptList(j).SwitchTimes, ...
            0, StartDaySkipTime, onlylongepoch);
        
        ffmeans_byday = structbird.All_FFdaymean;
        days = structbird.All_Daynum;
        
    end
    
    assert(all(days == alldaynums(allexptnum==j)), 'stimeffects and day means are not aligned ...');
    % ================= get slope of change in ff over day
    [~,~,~,~,~,SummaryStats]=lt_regress(ffmeans_byday,days, 0);
    
    SlopesAll_ffvsday = [SlopesAll_ffvsday; [SummaryStats.slope SummaryStats.slopeCI]];
    
    % ================= GET mean stim effect (ff, stim minus no stim)
    stimeffect = mean(allffstim(allexptnum==j));
    a = lt_sem(allffstim(allexptnum==j));
    StimEffectMean = [StimEffectMean; [stimeffect stimeffect-a stimeffect+a]];
     
    
    % ======== get day by day
    DayStats_StimEffect_MeanFF_DayNum = [DayStats_StimEffect_MeanFF_DayNum; ...
        [allffstim(allexptnum==j) ffmeans_byday days]];
end


lt_figure; hold on;
title('association experiments (50% stim)');
xlabel('ff change over days (slope)');
ylabel('mean [over days] effect of stim (ff)');
lt_plot(SlopesAll_ffvsday(:,1), StimEffectMean(:,1));
% --- plot CI
for j=1:size(SlopesAll_ffvsday,1)
   x = SlopesAll_ffvsday(j, :);
    y = StimEffectMean(j,:);
    
    line([x(2) x(3)], [y(1) y(1)], 'Color', 'k');
    line([x(1) x(1)], [y(2) y(3)], 'Color', 'k');
    
    tmp = [ExptList(j).birdname '-' ExptList(j).exptID];
    lt_plot_text(x(1)+0.3, y(1)+0.3, tmp, 'r');
end
lt_plot_zeroline;
lt_plot_zeroline_vert;



% ============================== IS THERE STRONGER STIM EFFECT ON DAYS WITH
% STRONGER LAERNIG?
lt_figure; hold on;

for j=1:length(DayStats_StimEffect_MeanFF_DayNum)
   
    dat = DayStats_StimEffect_MeanFF_DayNum{j};
    
    meanff = nan(max(dat(:,3)),1);
    stimeffect = nan(max(dat(:,3)),1);
    
    meanff(dat(:,3)) = dat(:,2);
    stimeffect(dat(:,3)) = dat(:, 1);
    
    % === for each day, get a metric of local slope
    if (1) % surroung
        locallearn = meanff(3:end) - meanff(1:end-2); % day n+1 minus n-1;
        stimeffect = stimeffect(2:end-1);
    elseif (0) % stim predict change net day?
        locallearn = meanff(2:end) - meanff(1:end-1); % day n+1 minus n-1;
        stimeffect = stimeffect(1:end-1);
    else % stim "amplifying" recent learning?
        locallearn = meanff(2:end) - meanff(1:end-1); % day n+1 minus n-1;
        stimeffect = stimeffect(2:end);
    end
	lt_subplot(4,2,j); hold on;
    title([ExptList(j).birdname '-' ExptList(j).exptID]);
    xlabel('stim effect(each day');
    ylabel('local learn (each day, n+1 - n-1)');
%     plot(stimeffect, locallearn, 'ok');
    lt_regress(locallearn, stimeffect, 1)
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end
























