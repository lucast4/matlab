function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S7(DATBYREND, TrialStruct, ...
    birdtoplot, expttotplot, mintime_fromref, maxtime_fromref, xedges)
%% lt 6/13/18 - plots for individual experiment (all syls)
dotrain=1;

%%

if expttotplot>length(TrialStruct.birds(birdtoplot).exptnum)
    disp('NO - this exptnum is too large, doesnt exist');
    return
end


indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot;
if ~any(indsthis)
    disp('NO DATA - likely is LMAN and not extracted ..')
    return
end


%%
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

xcenters = xedges(1:end-1)+diff(xedges)/2;


birdname = TrialStruct.birds(birdtoplot).birdname;
exptname = TrialStruct.birds(birdtoplot).exptnum(expttotplot).exptname;
numsyls = length(TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum);
for ss=1:numsyls
    syltoplot = ss;
    sylname = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(syltoplot).syl;
    
    istarg = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(syltoplot).INFO_istarget;
    issame = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(syltoplot).INFO_similar;
    if istarg==1 | issame==0
        continue
    end
    
    % ############################ 1) DENSITY IN REF PERIOD
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);        xlabel('N nontarg in ref');
    ylabel('N targ in ref');
    title('num rends in ref period');
    indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot & ...
        DATBYREND.Sylnum==syltoplot & DATBYREND.IsDurTrain==dotrain;
    
    nTARG = DATBYREND.Density_targ(indsthis);
    nNONTARG = DATBYREND.Density_nontarg(indsthis);
    hiDens = DATBYREND.Density_isHigh(indsthis);
    
    % --- add random jitter before plot
    nTARG = nTARG -0.3 + 0.6*rand(size(nTARG));
    nNONTARG = nNONTARG -0.3 + 0.6*rand(size(nNONTARG));
    % --- hi dens
    hd = 1;
    pcol = 'r';
    plot(nNONTARG(hiDens==hd), nTARG(hiDens==hd), 'x', 'Color', pcol);
    
    % --- lo dens
    hd = 0;
    pcol = 'b';
    plot(nNONTARG(hiDens==hd), nTARG(hiDens==hd), 'x', 'Color', pcol);
    
    
    % ############################ 2) ALL RENDITIONS, OVERLAY STUFF
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '-' exptname '-' sylname]);
    % 1) PLOT ALL TRIALS
    % - targ
    indtargsyl = find([TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum.INFO_istarget]);
    tvals_targ = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(indtargsyl).Tvals;
    ffvals_targ = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(indtargsyl).FFvals;
    plot(tvals_targ, ffvals_targ, '.k');
    
    % - nontarg
    indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot & ...
        DATBYREND.Sylnum==syltoplot;
    
    hiDens_tmp = DATBYREND.Density_isHigh(indsthis);
    timedevNT_tmp = DATBYREND.Time_dev(indsthis);
    tvals_tmp = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(syltoplot).Tvals;
    ffvals_tmp = TrialStruct.birds(birdtoplot).exptnum(expttotplot).sylnum(syltoplot).FFvals;
    assert(size(hiDens_tmp,1) == size(tvals_tmp,1), 'asfasd');
    
    plot(tvals_tmp, ffvals_tmp, '.b'); % ALL RENDS
    plot(tvals_tmp(hiDens_tmp==1), ffvals_tmp(hiDens_tmp==1), '.r'); % HI DENS
    plot(tvals_tmp(hiDens_tmp==0), ffvals_tmp(hiDens_tmp==0), '.c'); % HI DENS
    
    wnon = TrialStruct.birds(birdtoplot).exptnum(expttotplot).WNontime;
    line([wnon wnon], ylim, 'Color', 'k');
    
    % 2) FOR EACH RENDITION, PLOT i) ref window, ii) window with data'
    for kk=1:length(tvals_tmp)
        
        if isempty(timedevNT_tmp{kk})
            continue
        end
        
        tthistmp = tvals_tmp(kk);
        ffthistmp = ffvals_tmp(kk);
        
        % ----- reference window
        if (0)
            flankright = tthistmp+flanktime_targ/(60*24);
            flankleft = tthistmp-flanktime_targ/(60*24);
            line([flankleft flankright], [ffthistmp ffthistmp], ...
                'Color', 'c')
        end
        
        % ---- WINDOW for collecting all values precedinga dn
        % post
        if (0)
            datedges = tthistmp + twind_plot./(24);
            line([datedges], [ffthistmp ffthistmp], ...
                'Color', [0.6 0.6 0.6], 'LineStyle', '--');
        end
        
        % ----- window for counting number of rends [in plot
        % that will make in next section]
        %         edgeright = tthistmp + edge_time/(60*24);
        %         line([tthistmp edgeright], [ffthistmp ffthistmp], ...
        %             'Color', 'm', 'LineStyle', ':')
        
        % ---- line for window in which collect data
        winddat = [tthistmp+mintime_fromref/(60*24) tthistmp+maxtime_fromref/(60*24)];
        line(winddat, [ffthistmp ffthistmp], 'Color', [0.7 0.7 0.7]);
        line([tthistmp winddat(2)], [ffthistmp ffthistmp], 'Color', [0.7 0.7 0.7], ...
            'LineStyle', ':');
        
        
        % ---- plot number of left and right rends collected
        if (0)
            nrighttmp = sum(timedevNT_tmp{kk} > 0);
            nmidtmp = sum(timedevNT_tmp{kk} == 0);
            nlefttmp = sum(timedevNT_tmp{kk} < 0);
            lt_plot_text(edgeright, ffthistmp, ...
                ['n=' num2str(nlefttmp) '-' num2str(nmidtmp) '-' num2str(nrighttmp)], [0.6 0.6 0.6], 5);
        end
    end
    
    
    
    % ############################ 2) ff dev, locked to ref [TRAIN]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('DUR TRAIN');
    ylabel('ffdev');
    xlabel('time from ref');
    nmin = 200; % datapoints
    
    % ---- HI DENS
    hidens = 1;
    pcol = 'r';
    % -
    indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot & ...
        DATBYREND.Sylnum==syltoplot & DATBYREND.Density_isHigh==hidens & DATBYREND.IsDurTrain==1;
    tdev = DATBYREND.Time_dev(indsthis);
    ffdev = DATBYREND.FF_dev(indsthis);
    tdev = cell2mat(tdev)*(24*60);
    ffdev = cell2mat(ffdev);
    if length(tdev)>nmin
        indsrand = randperm(length(tdev), nmin);
        tdev = tdev(indsrand);
        ffdev = ffdev(indsrand);
    end
    plot(tdev, ffdev, 'x', 'Color', pcol);
    % -- overlay means
    tdev_binned = discretize(tdev, xedges);
    [ffmean, ffsem] = grpstats(ffdev, tdev_binned, {'mean', 'sem'});
    xmean = xcenters(unique(tdev_binned));
    lt_plot(xmean, ffmean, {'Errors', ffsem, 'Color', 'k', 'MarkerFaceColor', pcol});
    
    % ---- LO DENS
    hidens = 0;
    pcol = 'b';
    % -
    indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot & ...
        DATBYREND.Sylnum==syltoplot & DATBYREND.Density_isHigh==hidens & DATBYREND.IsDurTrain==1;
    tdev = DATBYREND.Time_dev(indsthis);
    ffdev = DATBYREND.FF_dev(indsthis);
    tdev = cell2mat(tdev)*(24*60);
    ffdev = cell2mat(ffdev);
    if length(tdev)>nmin
        indsrand = randperm(length(tdev), nmin);
        tdev = tdev(indsrand);
        ffdev = ffdev(indsrand);
    end
    plot(tdev, ffdev, 'x', 'Color', pcol);
    % -- overlay means
    tdev_binned = discretize(tdev, xedges);
    [ffmean, ffsem] = grpstats(ffdev, tdev_binned, {'mean', 'sem'});
    xmean = xcenters(unique(tdev_binned));
    lt_plot(xmean, ffmean, {'Errors', ffsem, 'Color', 'k', 'MarkerFaceColor', pcol});
    
    % ---
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    
    % ############################ 2) ff dev, locked to ref [BASELINE]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('DUR BASE');
    ylabel('ffdev');
    xlabel('time from ref');
    nmin = 200; % datapoints
    
    % ---- HI DENS
    hidens = 1;
    pcol = 'r';
    % -
    indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot & ...
        DATBYREND.Sylnum==syltoplot & DATBYREND.Density_isHigh==hidens & DATBYREND.IsDurTrain==0;
    tdev = DATBYREND.Time_dev(indsthis);
    ffdev = DATBYREND.FF_dev(indsthis);
    tdev = cell2mat(tdev)*(24*60);
    ffdev = cell2mat(ffdev);
    if length(tdev)>nmin
        indsrand = randperm(length(tdev), nmin);
        tdev = tdev(indsrand);
        ffdev = ffdev(indsrand);
    end
    plot(tdev, ffdev, 'x', 'Color', pcol);
    % -- overlay means
    tdev_binned = discretize(tdev, xedges);
    [ffmean, ffsem] = grpstats(ffdev, tdev_binned, {'mean', 'sem'});
    xmean = xcenters(unique(tdev_binned));
    lt_plot(xmean, ffmean, {'Errors', ffsem, 'Color', 'k', 'MarkerFaceColor', pcol});
    
    % ---- LO DENS
    hidens = 0;
    pcol = 'b';
    % -
    indsthis = DATBYREND.Birdnum==birdtoplot & DATBYREND.Exptnum==expttotplot & ...
        DATBYREND.Sylnum==syltoplot & DATBYREND.Density_isHigh==hidens & DATBYREND.IsDurTrain==0;
    tdev = DATBYREND.Time_dev(indsthis);
    ffdev = DATBYREND.FF_dev(indsthis);
    tdev = cell2mat(tdev)*(24*60);
    ffdev = cell2mat(ffdev);
    if length(tdev)>nmin
        indsrand = randperm(length(tdev), nmin);
        tdev = tdev(indsrand);
        ffdev = ffdev(indsrand);
    end
    plot(tdev, ffdev, 'x', 'Color', pcol);
    % -- overlay means
    tdev_binned = discretize(tdev, xedges);
    [ffmean, ffsem] = grpstats(ffdev, tdev_binned, {'mean', 'sem'});
    xmean = xcenters(unique(tdev_binned));
    lt_plot(xmean, ffmean, {'Errors', ffsem, 'Color', 'k', 'MarkerFaceColor', pcol});
    
    % ---
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
end