function lt_neural_Coher_LFP_plottrials(SwitchCohStruct, OUTSTRUCT, SwitchStruct, ...
    PARAMS, birdtoplot, expttoplot, swnum, motiftoplot, chanpairtoplot, fs, ...
    Nplot, bregiontoget, plotonlyLFP)

%% lt 10/31/18 - plots individual trials LFP,coherence, sp[ectrogram.
% --- takes N example trials, pre and post learning, for a given motif

% birdtoplot = 'pu69wh78';
% expttoplot = 'RALMANlearn2';
% swnum = 1;
% motiftoplot = ''; % ------ WILL plot target syl, unless say otherwise.
% chanpairtoplot = []; % if empty, will pick first pari in list.

% Nplot = 5; % how any trials each for base and WN?


%%
% plotON=0;
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;
numbirds = length(SwitchCohStruct.bird);


figcount=1;
if plotonlyLFP==1
    subplotcols = 2;
    subplotrows=4;
else
    subplotcols=6;
    subplotrows=3;
end
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ===== chronux params

lt_switch_chronux(1);
movingwin = [0.1 0.01];

params = struct;
params.fpass = [1/movingwin(1) 150];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500; % hard coded fs for LFP;

lfpband_lo = 25;
lfpband_hi = 35;

ffhilim = 90; % for plots;


% ======== for power spectra
twind_spectrum = [-0.08 0]; % relative to syl onset

%%


% ======================
for i=1:numbirds
    bthis = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        foundatarg =0;
        
        for ss=1:numswitch
            
            if ~(strcmp(birdtoplot, bthis) & strcmp(expttoplot, ename) & ss==swnum)
                continue
            end
            
            nummotifs = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
            datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss);
            if isempty(datthis.motifnum)
                continue
            end
            
            % =========== for this switch, get edge times and switch time
            %             tstart = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            %             tswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            %             tend = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            
            targsyls = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{1:2:end}};
            targdirs = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}];
            motiflist = {datthis.motifnum.motifname};
            motiflist = motiflist(~cellfun(@isempty, motiflist));
            
            sylssame = lt_neural_QUICK_ExtractSameType(motiflist, targsyls);
            
            
            if isempty(motiftoplot)
                motiftoplot = targsyls{1};
            end
            
            mm = find(strcmp({datthis.motifnum.motifname}, motiftoplot));
            assert(length(mm)==1, 'motif not exist?');
            
            
            % =============================== DATA (SWITCH COH STRUCT)
            if isempty(chanpairtoplot)
                chanpairtoplot = datthis.motifnum(mm).chanpair(1,:);
                tmp = find(strcmp(datthis.motifnum(mm).bregionpair, bregiontoget));
                chanpairtoplot = datthis.motifnum(mm).chanpair(randperm(length(tmp),1),:);
            end
            
            % ============== COHERENCE DATA
            tmp2 = chanpairtoplot == datthis.motifnum(mm).chanpair;
            indchanpair = find(all(tmp2'));
            
            if isfield(datthis.motifnum(mm), 'cohmat')
                cohmat = squeeze(datthis.motifnum(mm).cohmat(:,:, :, indchanpair));
            else
                cohmat = lt_neural_LFP_QUICK_loadProcessDat(...
                    [datthis.motifnum(mm).fileprefix '/Coh' datthis.motifnum(mm).filesuffix], indchanpair);
            end
            %             tvals = datthis.motifnum(mm).tvals;
            bregionpair = datthis.motifnum(mm).bregionpair{indchanpair};
            t_coh = PARAMS.tbins;
            ff_coh = PARAMS.ffbins;
            
            % ============== LFP DATA
            indlfpchan = find(ismember(datthis.motifnum(mm).lfpall_chans, chanpairtoplot));
            
            lfpmat1 = cell2mat(datthis.motifnum(mm).lfpall(:,indlfpchan(1))')';
            lfpmat2 = cell2mat(datthis.motifnum(mm).lfpall(:,indlfpchan(2))')';
            
            lfpchan1 = datthis.motifnum(mm).lfpall_chans(indlfpchan(1));
            lfpchan2 = datthis.motifnum(mm).lfpall_chans(indlfpchan(2));
            
            t_lfp = datthis.motifnum(mm).t_lfp;
            
            % ================= high pass fuilter the lfp?
            %                 [filt_b,filt_a]=butter(2,[2*2/fs], 'high');
            
            
            %     % === filter
            %     datfilt=filtfilt(filt_b,filt_a,double(lfpmat1));
            %
            %             lfpmat11 = lt_neural_filter(double(lfpmat1), fs, 0, 2, 400);
            
            %% ============== METADATA
            % ---------- FROM OUTSTRUCT
            if isempty(OUTSTRUCT)
                disp('NOTE: will not separate by WN and base, since did not enter OUTSTRUCT..')
                pause;
                indsbase = 1:size(cohmat,3);
                indsWN = 1:size(cohmat,3);
            else
                
                indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss ...
                    & OUTSTRUCT.motifnum==mm);
                indsbase = OUTSTRUCT.indsbase_epoch{indsthis(1)};
                indsWN = OUTSTRUCT.indsWN_epoch{indsthis(1)};
            end
            
            
            % ----------- COLLECT MEAN COHERENCE PRE AND POST SWITCH
            %             cohmat = datthis.motifnum(mm).cohmat;
            %             chanpairs = datthis.motifnum(mm).chanpair;
            %
            %             % ------ take mean across all channels paiors
            %             if averagechanpairs ==1
            %                 cohmat = nanmean(cohmat, 4);
            %             end
            %
            % --- TARG? what type of syl is this?
            if any(strcmp(targsyls, motiftoplot))
                istarg =1;
                learndir = targdirs(find(strcmp(targsyls, motiftoplot)));
                titlecol = 'r';
            else
                istarg = 0;
                learndir = [];
                titlecol = 'k';
            end
            
            % --- SAME?
            issame = any(strcmp(sylssame, motiftoplot));
            
            predur_lfp = -t_lfp(1); % i.e. combined motif_predur + flank time
            
            
            %% ====== plots
            
            % ======== get random trials base and WN
            indsbase_shuff = indsbase(randperm(length(indsbase), Nplot));
            indsWN_shuff = indsWN(randperm(length(indsWN), Nplot));
            
            
            % ====== PLOT EACH TRIAL (LEFT: BASE... RIGHT: WN)
            for j=1:Nplot
                
                % ========================================== BASELINE
                indtrial = indsbase_shuff(j);
                plotname = 'BASE';
                
                lt_neural_Coher_LFP_plottrials_sub1;
                
                
                % ========================================== WN
                indtrial = indsWN_shuff(j);
                plotname = 'WN';
                
                lt_neural_Coher_LFP_plottrials_sub1;
                
            end
            
            %% ===================== SUMMARY PLOT (BASE AND WN)
            lt_figure; hold on;
            
            % ================================= MEAN LFP
            % ----------------- CHANNEL 1
            lt_subplot(4,2,1); hold on;
            title(['chan1' bregionpair(1:4)]);
            
            lfpmat_filt = lt_neural_filter(double(lfpmat1'), fs, 0, lfpband_lo, lfpband_hi)';
            
            % - base
            lfpmean = mean(lfpmat_filt(indsbase,:),1);
            lfpsem = lt_sem(lfpmat_filt(indsbase,:));
            shadedErrorBar(t_lfp, lfpmean, lfpsem, {'Color', 'k'},1);
            
            % - WN
            lfpmean = mean(lfpmat_filt(indsWN,:),1);
            lfpsem = lt_sem(lfpmat_filt(indsWN,:));
            shadedErrorBar(t_lfp, lfpmean, lfpsem, {'Color', 'm'},1);
            
            
            % ----------------- CHANNEL 2
            lt_subplot(4,2,2); hold on;
            title(['chan2' bregionpair(end-3:end)]);
            
            lfpmat_filt = lt_neural_filter(double(lfpmat2'), fs, 0, lfpband_lo, lfpband_hi)';
            
            % - base
            lfpmean = mean(lfpmat_filt(indsbase,:),1);
            lfpsem = lt_sem(lfpmat_filt(indsbase,:));
            shadedErrorBar(t_lfp, lfpmean, lfpsem, {'Color', 'k'},1);
            
            % - WN
            lfpmean = mean(lfpmat_filt(indsWN,:),1);
            lfpsem = lt_sem(lfpmat_filt(indsWN,:));
            shadedErrorBar(t_lfp, lfpmean, lfpsem, {'Color', 'm'},1);
            
            % ===================================== MEAN SPECTROGRAM
            if plotonlyLFP==0
                paramstmp = params;
                paramstmp.trialave = 1;
                
                % -------- BASE (CHAN1)
                lt_subplot(4,2,3); hold on;
                title('BASE (chan1)');
                [S, t, f] = mtspecgramc(lfpmat1(indsbase,:)', movingwin, paramstmp);
                t = t-predur_lfp;
                imagesc(t, f, 10*log10(S)');
                colorbar
                axis tight;
                ylim([0 ffhilim]);
                
                % -------- BASE (CHAN2)
                lt_subplot(4,2,4); hold on;
                title('BASE (chan2)');
                [S, t, f] = mtspecgramc(lfpmat2(indsbase,:)', movingwin, paramstmp);
                t = t-predur_lfp;
                imagesc(t, f, 10*log10(S)');
                colorbar
                axis tight;
                ylim([0 ffhilim]);
                
                % -------- WN (CHAN1)
                lt_subplot(4,2,5); hold on;
                title('WN (chan1)');
                [S, t, f] = mtspecgramc(lfpmat1(indsWN,:)', movingwin, paramstmp);
                t = t-predur_lfp;
                imagesc(t, f, 10*log10(S)');
                colorbar
                axis tight;
                ylim([0 ffhilim]);
                
                
                % -------- WN (CHAN2)
                lt_subplot(4,2,6); hold on;
                title('WN (chan2)');
                [S, t, f] = mtspecgramc(lfpmat2(indsWN,:)', movingwin, paramstmp);
                t = t-predur_lfp;
                imagesc(t, f, 10*log10(S)');
                colorbar
                axis tight;
                ylim([0 ffhilim]);
                
                
                % NOTE: This first normalizes each trial so that power peak is
                % at 1, and then averages.
                %             % ========================= SPECTRA
                %             paramstmp = params;
                %             paramstmp.err = [2 0.05];
                %             paramstmp.trialave = 0;
                %
                %             % ------- CHAN1 (base and WN)
                %             lt_subplot(4,2,7); hold on;
                %             title('chan1 (base=k, wn=r)');
                %             lfpmatthis = lfpmat1;
                %
                %             % base
                %             [S, f, Serr] = mtspectrumc(lfpmatthis(indsbase,:)', paramstmp);
                %             Snorm = S./repmat(max(S, [], 1), size(S,1), 1);
                %             %             plot(f, Snorm, '-k');
                %             shadedErrorBar(f, mean(Snorm,2), lt_sem(Snorm'), {'Color', 'k'},1);
                %
                %             % WN
                %             [S, f, Serr] = mtspectrumc(lfpmatthis(indsWN,:)', paramstmp);
                %             Snorm = S./repmat(max(S, [], 1), size(S,1), 1);
                %             shadedErrorBar(f, mean(Snorm,2), lt_sem(Snorm'), {'Color', 'r'},1);
                %
                %             % ------- CHAN2 (base and WN)
                %             lt_subplot(4,2,8); hold on;
                %             title('chan2 (base=k, wn=r)');
                %             lfpmatthis = lfpmat2;
                %
                %             % base
                %             [S, f, Serr] = mtspectrumc(lfpmatthis(indsbase,:)', paramstmp);
                %             Snorm = S./repmat(max(S, [], 1), size(S,1), 1);
                %             %             plot(f, Snorm, '-k');
                %             shadedErrorBar(f, mean(Snorm,2), lt_sem(Snorm'), {'Color', 'k'},1);
                %
                %             % WN
                %             [S, f, Serr] = mtspectrumc(lfpmatthis(indsWN,:)', paramstmp);
                %             Snorm = S./repmat(max(S, [], 1), size(S,1), 1);
                %             shadedErrorBar(f, mean(Snorm,2), lt_sem(Snorm'), {'Color', 'r'},1);
                %
                %             if (0)
                %                figure; hold on;
                %                plot(f, S, '-k');
                %                figure; hold on;
                %                lt_plot_histogram(10*log10(S(5,:)))
                %             end
                
                % ========================= SPECTRA
                paramstmp = params;
                paramstmp.err = [2 0.05];
                paramstmp.trialave = 1;
                
                % ------- CHAN1 (base and WN)
                lt_subplot(4,2,7); hold on;
                title('chan1 (base=k, wn=r)');
                lfpmatthis = lfpmat1;
                
                % base
                [S, f, Serr] = mtspectrumc(lfpmatthis(indsbase,:)', paramstmp);
                plot_vector(S, f, 'l', Serr, 'k');
                % WN
                [S, f, Serr] = mtspectrumc(lfpmatthis(indsWN,:)', paramstmp);
                plot_vector(S, f, 'l', Serr, 'r');
                
                % ------- CHAN2 (base and WN)
                lt_subplot(4,2,8); hold on;
                title('chan2 (base=k, wn=r)');
                lfpmatthis = lfpmat2;
                
                % base
                [S, f, Serr] = mtspectrumc(lfpmatthis(indsbase,:)', paramstmp);
                plot_vector(S, f, 'l', Serr, 'k');
                % WN
                [S, f, Serr] = mtspectrumc(lfpmatthis(indsWN,:)', paramstmp);
                plot_vector(S, f, 'l', Serr, 'r');
            end
        end
    end
end

lt_switch_chronux(0);
