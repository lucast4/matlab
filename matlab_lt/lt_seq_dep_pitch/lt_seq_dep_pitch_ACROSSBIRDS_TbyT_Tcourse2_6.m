function DATTMP = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
    TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn)

%%
% istrain = 1;
% iscatch = 0;
% 
% istrain = [0 1];
% iscatch = 0;
% 
% istrain = 0;
% iscatch = [0 1];
% ------ 


DATTMP = struct;

%%
for i=1:length(Inds_sylcounter)
   
    ss = Inds_sylcounter(i);
    
    indsthis = DATBYREND.Sylcounter==ss & ismember(DATBYREND.IsDurTrain, istrain) ...
    & ismember(DATBYREND.IsCatch, iscatch) ...        
    & ~cellfun(@isempty, DATBYREND.FF_dev);
    
    if ~any(indsthis)
        continue
    end
    
    ffthis = DATBYREND.FF_dev(indsthis);
    tthis = DATBYREND.Time_dev(indsthis);
    if use_targ_locallearn==1
    learnlocal = DATBYREND.LearnLocal_targ(indsthis);
    else
        learnlocal = DATBYREND.LearnLocal(indsthis);
    end
    bname = TrialStruct.birds(unique(DATBYREND.Birdnum(indsthis))).birdname;
    ename = TrialStruct.birds(unique(DATBYREND.Birdnum(indsthis))).exptnum(unique(DATBYREND.Exptnum(indsthis))).exptname;
    
    % ############################## COLLECT STUFF
    % ========= 1) bin all data by magnitude of local learnig
    
    
    % ========= 2) collect [locallearn, first trial dev]
    functmp = @(x)(x(1));
    ffdev_first = cellfun(functmp, ffthis);
    tdev_first = cellfun(functmp, tthis).*(24*60);
    assert(all(cellfun(functmp, tthis))>0, 'i cant use this func to collect first deviation of pitch');
    DATTMP(i).learnlocal = learnlocal;
    DATTMP(i).ffdev_first = ffdev_first;
    DATTMP(i).tdev_first = tdev_first;
    
    
    if (0)
    % ======================== SEPARATE INTO HIGH AND LOW LOCAL LAERNING
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    %         title(['targ']);
    title([bname '-' ename]);
    xlabel('local learn (e.. trial n+1)');
    ylabel('next trial learn (e.g. trial n+2)');
%     plot(learnlocal, 
    
    
    
    % ======================== SEPARATE INTO HIGH AND LOW LOCAL LAERNING
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    %         title(['targ']);
    title([bname '-' ename]);
    xlabel('time dev');
    ylabel('ff dev');
    % - high
    pcol = 'r';
    indtmp = learnlocal>0;
    ttmp = cell2mat(tthis(indtmp))*60*24;
    fftmp = cell2mat(ffthis(indtmp));
    plot(ttmp, fftmp, 'x', 'Color', pcol);

    % - lo
    pcol = 'b';
    indtmp = learnlocal<0;
    ttmp = cell2mat(tthis(indtmp))*60*24;
    fftmp = cell2mat(ffthis(indtmp));
    plot(ttmp, fftmp, 'x', 'Color', pcol);
    end
    
end

% if length(DATTMP)<length(Inds_sylcounter)
%    DATTMP( 
% end
