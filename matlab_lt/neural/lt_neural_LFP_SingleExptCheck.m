%% ++++++++++++++ 1) COHSCALAR, OVER TRIALS
close all;

for i=1:length(SummaryStruct.birds)
    %     birdplot = 'pu69wh78';
    birdplot = SummaryStruct.birds(i).birdname;
    for ii=1:length(SwitchStruct.bird(i).exptnum)
        exptplot = SummaryStruct.birds(i).exptnum_pop(ii).exptname;
        % exptplot = 'RALMANlearn1';
        for ss=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist)
            swplot = ss;
            % swplot = 1;
            
            motifplot = []; % [string] leave blank for target
            
            %             i = find(strcmp({SummaryStruct.birds.birdname}, birdplot));
            %             ii = find(strcmp({SummaryStruct.birds(i).exptnum_pop.exptname}, exptplot));
            %             ss = swplot;
            
            % ===========
            
            % --- inds for this experiment [target syl]
            indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
                OUTSTRUCT.istarg==1);
            if ~any(indsthis)
                continue
            end
            
            % --- target syl motif number?
            mm = unique(OUTSTRUCT.motifnum(indsthis));
            assert(length(mm)==1, 'multipel targets?');
            
            
            % ######################3 EXTRACT DATA AND PLOT
            chanpairs = OUTSTRUCT.chanpair(indsthis, :);
            
            tvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).tvals;
            ffvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).ffvals;
            cohscal_allpairs = OUTSTRUCT.cohscal(indsthis);
            indsbase = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch;
            indsWN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch;
            indsbase_all = OUTSTRUCT.indsbase{indsthis(1)};
            indsWN_all = OUTSTRUCT.indsWN{indsthis(1)};
            
            
            % ================ PLOT FF
            figcount=1;
            subplotrows=4;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots; hsplot];
            title([birdplot '-' exptplot '-sw' num2str(swplot)]);
            ylabel('ff');
            plot(tvals, ffvals, 'ok');
            axis tight;
            
            % ------ mark baseline
            line([min(tvals(indsbase_all)) min(tvals(indsbase_all))], ylim, 'Color', 'b');
            line([max(tvals(indsbase_all)) max(tvals(indsbase_all))], ylim, 'Color', 'b');
            YLIM = ylim;
            lt_neural_QUICK_PlotSylPatches(tvals(indsbase(1)), tvals(indsbase(end)), YLIM, 0, 'b');
            % ------ mark WN end
            line([min(tvals(indsWN_all)) min(tvals(indsWN_all))], ylim, 'Color', 'r');
            line([max(tvals(indsWN_all)) max(tvals(indsWN_all))], ylim, 'Color', 'r');
            YLIM = ylim;
            lt_neural_QUICK_PlotSylPatches(tvals(indsWN(1)), tvals(indsWN(end)), YLIM, 0, 'r');
            
            
            % ================= PLOT coherence
            for cc=1:size(chanpairs,1)
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots; hsplot];
                title(['chan ' num2str(chanpairs(cc,:))]);
                plot(tvals, cohscal_allpairs{cc}, 'ok');
                axis tight;
                ylim([0 1]);
                
                % ------ mark baseline
                line([min(tvals(indsbase_all)) min(tvals(indsbase_all))], ylim, 'Color', 'b');
                line([max(tvals(indsbase_all)) max(tvals(indsbase_all))], ylim, 'Color', 'b');
                YLIM = ylim;
                lt_neural_QUICK_PlotSylPatches(tvals(indsbase(1)), tvals(indsbase(end)), YLIM, 0, 'b');
                % ------ mark WN end
                line([min(tvals(indsWN_all)) min(tvals(indsWN_all))], ylim, 'Color', 'r');
                line([max(tvals(indsWN_all)) max(tvals(indsWN_all))], ylim, 'Color', 'r');
                YLIM = ylim;
                lt_neural_QUICK_PlotSylPatches(tvals(indsWN(1)), tvals(indsWN(end)), YLIM, 0, 'r');
            end
            
            
            % ================= PLOT CHANGE IN COHERENCE ACROSS ALL CHANNELS
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('means within shaded periods');
            ylabel('coherence');
            xlabel('pre -- post');
            for cc=1:size(chanpairs,1)
                cohmean = [mean(cohscal_allpairs{cc}(indsbase)) mean(cohscal_allpairs{cc}(indsWN))];
                plot([1 2], cohmean, '-ok');
                lt_plot_text(2.2, cohmean(2), num2str(chanpairs(cc,:)));
            end
            xlim([0 3]);
            
            
            % ================== get xcorr between ff and coherence
            % scalars.
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            %             title('means within shaded periods');
            %             ylabel('coherence');
            %             xlabel('pre -- post');
            xlabel('coh leads -- ff leads [trials]');
            ylabel('xcov');
            for cc=1:size(chanpairs,1)
                cohthis = cohscal_allpairs{cc}(indsbase(1):indsWN(end));
                ffthis = ffvals(indsbase(1):indsWN(end));
                [ccov, lags] = xcov(cohthis', ffvals');
                plot(lags, ccov, 'k');
            end
            xlim([-20 20]);
            lt_plot_zeroline;
            
            % ------------------
            linkaxes(hsplots, 'x');
        end
    end
end


%% ++++++++++++++ 2) RAW NEURAL DATA
%         birdplot = 'pu69wh78';
%         exptplot = 'RALMANlearn1';
birdplot = 'wh44wh39';
exptplot = 'RALMANlearn3';
swplot = 1;
motifplot = []; % [string] leave blank for target

i = find(strcmp({SummaryStruct.birds.birdname}, birdplot));
ii = find(strcmp({SummaryStruct.birds(i).exptnum_pop.exptname}, exptplot));
ss = swplot;

indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
    OUTSTRUCT.istarg==1);

% --- target syl motif number?
mm = unique(OUTSTRUCT.motifnum(indsthis));
assert(length(mm)==1, 'multipel targets?');


% ######################3 EXTRACT DATA AND PLOT
chanpairs = OUTSTRUCT.chanpair(indsthis, :);

tvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).tvals;
ffvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).ffvals;
cohscal_allpairs = OUTSTRUCT.cohscal(indsthis);
indsbase = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch;
indsWN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch;
indsbase_all = OUTSTRUCT.indsbase{indsthis(1)};
indsWN_all = OUTSTRUCT.indsWN{indsthis(1)};



% ======== 1) EXTRACT RAW DATA (NOT JUST LFP)
%     close all;

%     birdplot = 'pu69wh78';
%     exptplot = 'RALMANlearn1';
%     swplot = 1;
%     motifplot = []; % [string] leave blank for target
extrapad = 0.05; % seconds, pre and post...
[DatAll, t_onoff, fs, bregionlist, chanlist_toget, i, ii, mm] = ...
    lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, ...
    SwitchCohStruct, birdplot, exptplot, swplot, motifplot, ...
    extrapad);


% =========== GET LIST OF HIGH AND LOW COHERENCE TRIALS
%     tscalar = -0.045;
%     fscalar = 30;
%
%     a = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).fileprefix;
%     b = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).filesuffix;
%     cohdat = load([a '/Coh' b]);
%     cohdat = lt_neural_Coher_Cell2Mat(cohdat.CohAllTrials);
%     [~, ind_f] = min(abs(PARAMS.ffbins-fscalar));
%     [~, ind_t] = min(abs(PARAMS.tbins-tscalar));
%     cohscaltmp = squeeze(cohdat(ind_t, ind_f, :));

% ============ 2) FOR BASE AND WN, PLOT N TRIALS OF LFP, NEURAL
%     close all;
ntoplot = 10; % trials pre and post
filt_low = 25;
filt_hi = 35;
% filt_fs = fs;

savedir = ['/bluejay5/lucas/analyses/neural/LFP/FIGS_PlotEgRaw/' PARAMS.savemarker];
saveON = 0;

plotCohScalExtremes =1; % if 1, then plots example higha and low coherence tirals.
if ~exist('cohscaltmp', 'var')
    cohscaltmp = [];
end
lt_neural_LFP_PlotEgRaw(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
    i, ii, swplot, mm, ntoplot, filt_low, filt_hi, SwitchCohStruct, ...
    SwitchStruct, savedir, saveON, plotCohScalExtremes, cohscaltmp, ...
    chanpairs, cohscal_allpairs);


%% ======= 4) FOR A GIVEN TRIAL, WALK THROUGH STEPS OF CALCULATING COHERENCE
% NOTE: USES DATA FROM ABOVE.

% ++++++++++++++++++++++++ PARAMS
trialnum = 1;
chanpair = [14 17];
% --- FOR SCALAR EXTRACTION  - just sanity check, comapre to prvious
% extraction.
twind = [-0.09 -0.02];
fwind = [15 40];



% === extract LFP
setnum = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).neursetused;

lfpdat = LFPSTRUCT.bird(i).experiment(ii).setnum(setnum).motif(mm);

ind1 = lfpdat.Chanlist==chanpair(1);
ind2 = lfpdat.Chanlist==chanpair(2);

lfp1 = lfpdat.LFP_chanbytrial{ind1, trialnum};
lfp2 = lfpdat.LFP_chanbytrial{ind2, trialnum};

dat1 = DatAll{chanlist_toget==chanpair(1)}(:,trialnum);
dat2 = DatAll{chanlist_toget==chanpair(2)}(:,trialnum);

if (0)
    lt_figure; hold on;
    
    subplot(2,1,1); hold on;
    plot(lfp1, 'b');
    plot(lfp2, 'r');
    
    subplot(2,1,2); hold on;
    plot(dat1, 'b');
    plot(dat2, 'r');
    
end

% ==== COHERENCE USING THIS LFP
movingwin = [0.1 0.01];
params = struct;
params.fpass = [1/movingwin(1) 150];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500; % hard coded fs for LFP;
% params.Fs = 30000; % hard coded fs for LFP;

lt_switch_chronux(1);
[C,phi,S12,S1,S2,t,ffbins] = cohgramc(lfp1, lfp2, movingwin, params);
% [C,phi,S12,S1,S2,t,ffbins] = cohgramc(dat1, dat2, movingwin, params);
lt_switch_chronux(0);
t = t - PARAMS.motif_predur - extrapad;

lt_figure; hold on;
lt_neural_Coher_Plot(C, t, ffbins, 1, '', [0 1]);
axis tight;



% #####################  EXTRACT JUST COHERENCE SPECTRUM - do this with different
% length data
datlengthlist = [0.3 0.2 0.1 0.05]; % in sec
triallist = 1:50;
C_all = cell(length(triallist), length(datlengthlist));
S12_all = cell(length(triallist), length(datlengthlist));
S1_all = cell(length(triallist), length(datlengthlist));
S2_all = cell(length(triallist), length(datlengthlist));
Phi_all = cell(length(triallist), length(datlengthlist));
f_all =cell(length(triallist), length(datlengthlist));
for tr = 1:length(triallist)
    trialthis = triallist(tr);
    lfp1 = lfpdat.LFP_chanbytrial{ind1, trialthis};
    lfp2 = lfpdat.LFP_chanbytrial{ind2, trialthis};
    
    for kkk=1:length(datlengthlist)
        disp([trialthis kkk])
        datlength = datlengthlist(kkk);
        
        % --- prune data to get this long.
        
        lfp1_tmp = lfp1(1:datlength*1500);
        lfp2_tmp = lfp2(1:datlength*1500);
        
        lt_switch_chronux(1);
        params = struct;
        params.fpass = [1/movingwin(1) 150];
        w = 15; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
        tw = movingwin(1)*w;
        params.tapers = [tw 2*tw-1];
        params.Fs = 1500; % hard coded fs for LFP;
        
        
        [C,phi,S12,S1,S2,f] = coherencyc(lfp1_tmp, lfp2_tmp, params);
        lt_switch_chronux(0);
        
        
        
        % === collect C
        C_all{tr, kkk} = C';
        f_all{tr, kkk} = f;
        S12_all{tr, kkk} = S12;
        S1_all{tr, kkk} = S1;
        S2_all{tr, kkk} = S2;
    end
end

% ========== RECALCULATE COHERENCE BY HAND
Cohere_matlab = cell(1, length(datlengthlist));
for kkk=1:length(datlengthlist)
        datlength = datlengthlist(kkk);
        
    lfp1 = lfpdat.LFP_chanbytrial(ind1, triallist);
    % -- prune to correct length
    lfp1 = cellfun(@(x)x(1:datlength*1500), lfp1, 'UniformOutput', 0);
    lfp1 = cell2mat(cellfun(@transpose, lfp1, 'UniformOutput', 0));
    
    lfp2 = lfpdat.LFP_chanbytrial(ind2, triallist);
    % -- prune to correct length
    lfp2 = cellfun(@(x)x(1:datlength*1500), lfp2, 'UniformOutput', 0);
    lfp2 = cell2mat(cellfun(@transpose, lfp2, 'UniformOutput', 0));
    
    window = length(lfpdat.LFP_chanbytrial{ind1, triallist(1)});
    [c, f] = mscohere(lfp1, lfp2, window, 0.5*window, f_all{1, kkk}, 1500);
    
    Cohere_matlab{kkk} = c;
end


% ========== RECALCULATE COHERENCE BY HAND, across all trials multitaper
% method - i.e. N = trial x ntapers
% NOTE: taking mean(trials) of mean(tapers) is identical to taking
% mean(trials x tapers) - even for complex numbers.
Cohere_multitaper = cell(1, length(datlengthlist));
for kkk=1:length(datlengthlist)
        datlength = datlengthlist(kkk);
        
        lfp1 = lfpdat.LFP_chanbytrial(ind1, triallist);
        % -- prune to correct length
        lfp1 = cellfun(@(x)x(1:datlength*1500), lfp1, 'UniformOutput', 0);
        lfp1 = cell2mat(lfp1);
        
        lfp2 = lfpdat.LFP_chanbytrial(ind2, triallist);
        % -- prune to correct length
        lfp2 = cellfun(@(x)x(1:datlength*1500), lfp2, 'UniformOutput', 0);
        lfp2 = cell2mat(lfp2);
        
        params = struct;
        params.fpass = [1/movingwin(1) 150];
        w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
        tw = movingwin(1)*w;
        params.tapers = [tw 2*tw-1];
        params.Fs = 1500; % hard coded fs for LFP;
        params.trialave = 1;

        lt_switch_chronux(1);
        [C,phi,S12,S1,S2,f] = coherencyc(lfp1, lfp2, params);
        lt_switch_chronux(0);

                
    Cohere_multitaper{kkk} = C;
end

% ============= CALCULATE COHERENCE USING EXTRQACT SPECTRA
C_xtrial = cell(1,length(datlengthlist));
for dd=1:size(S12_all,2)
    s12this = cell2mat(cellfun(@transpose, S12_all(:,dd), 'UniformOutput', 0));
    s1this = cell2mat(cellfun(@transpose, S1_all(:,dd), 'UniformOutput', 0));
    s2this = cell2mat(cellfun(@transpose, S2_all(:,dd), 'UniformOutput', 0));
    
    s12this = mean(s12this,1)';
    s1this = mean(s1this,1)';
    s2this = mean(s2this,1)';
    c12 = abs(s12this./sqrt(s1this.*s2this));
    
    C_xtrial{dd} = c12;
end


% === PLOT
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for dd=1:size(C_all,2) % iterate over all data sizes
    Cthis = cell2mat(C_all(:,dd)); % get all trials for this datasize
    fthis = f_all{1,dd};
    cmean = mean(Cthis,1);
    csem = lt_sem(Cthis);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['data length ' num2str(datlengthlist(dd)) 'sec']);
    ylabel('b(x trial), g(x trial/tapers), r(matlab)');
    plot(fthis, Cthis', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(fthis, cmean, csem, {'Color', 'r'}, 1);
    
    % -- overlay aross trial coherence
    plot(fthis, C_xtrial{dd}, '--b', 'LineWidth', 2);
    axis tight
    
    % ---overlay matlab version
    plot(fthis, Cohere_matlab{dd}, '-r');
    
    % -- overlay multitaper over tapers/trials
    plot(fthis, Cohere_multitaper{dd}, 'g');
end
linkaxes(hsplots,'xy');




% ############################## EXTRACT SCALAR TO COMPARE TO PRVIOUS DATA
mean(mean(C(t>twind(1) & t<twind(2), ffbins>fwind(1) & ffbins<fwind(2))));

pairthis = all(chanpairs'==chanpair');
cohscal_allpairs{pairthis}(trialnum)



%% ======= FOR A GIVEN TRIAL, SHOW ALL SEPARATE TAPERS 
% 1) Each taper show raw data, show filtered data (30hz)
% 2) Show untapered data, filtered
% 3) Show coherence scalar.

% NOTE: first runt he script avbove, plotting individual trials...
trialthis = 869;
pairtoplot = [15 17];
twindtoget = [-0.1 0]; % relative syl onset
ftoplot = 30;
assert(ftoplot==30, 'need to adjust neural filter...');

% ==== what is coernee scalar?
chpairtmp = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).chanpair;
cohscal = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).cohscalar(trialthis, ...
    all(chpairtmp' == pairtoplot'));


% ==== extract lfp
setnum = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).neursetused;
lfpdat = LFPSTRUCT.bird(i).experiment(ii).setnum(setnum).motif(mm);

ind1 = lfpdat.Chanlist==pairtoplot(1);
ind2 = lfpdat.Chanlist==pairtoplot(2);

lfp1 = lfpdat.LFP_chanbytrial{ind1, trialthis};
lfp2 = lfpdat.LFP_chanbytrial{ind2, trialthis};

t = lfpdat.t_relons;

indsthis = t>twindtoget(1) & t<twindtoget(2);

% ==== segmented LFP DAT:
lfp1 = lfp1(indsthis);
lfp2 = lfp2(indsthis);
t = t(indsthis);

% ==== coherency params, used to extract appropriate tapered raw data.

movingwin = [0.1 0.01];
params = struct;
params.fpass = [1/movingwin(1) 150];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500; % hard coded fs for LFP;

[data1, data2, data_proj1, data_proj2, tapers, J12, f] = ...
    lt_neural_Coher_taperdata(lfp1,lfp2,params);


% ============================ PLOT
lt_figure; hold on;

% === get cross spectra for each taper
[~, indf] = min(abs(f-ftoplot));

nrow = 6;

% === raw data
lt_subplot(nrow, 3, 1); hold on;
title(['b(ch1), r(ch2), COH=' num2str(cohscal)]);
plot(t, lfp1, 'b');
plot(t, lfp2,  'r');

% === raw data (filtered)
lt_subplot(nrow, 3, 2); hold on;
lfp1_filt = lt_neural_filter(double(lfp1), 1500, 0, 28, 32);
lfp2_filt = lt_neural_filter(double(lfp2), 1500, 0, 28, 32);
plot(t, lfp1_filt, 'b');
plot(t, lfp2_filt,  'r');

% ==== plot cross spectru,
j = mean(J12(indf,:), 2);
lt_subplot(nrow, 3, 3); hold on;
plot(real(j), imag(j), 'ok');
xlim([-30 30]);
ylim([-30 30]);
lt_plot_zeroline;
lt_plot_zeroline_vert;




% ==== each taper
ntapers = size(data1,2);
for nt=1:ntapers
    
    % ---- raw data
    lt_subplot(nrow, 3, (3*nt-2)+3); hold on;
    title('b(ch1), r(ch2)');
    ylabel('raw, after taper');
    plot(t, data_proj1(:,nt), 'b');
    plot(t, data_proj2(:,nt),  'r');
    
    % -- overlay tapers
    plot(t, 50*tapers(:,nt), '-k');
    
    % ---- raw data (filtered)
    lt_subplot(nrow, 3, (3*nt-1)+3); hold on;
    plot(t, lt_neural_filter(double(data_proj1(:,nt)), 1500, 0, 28, 32), 'b');
    plot(t, lt_neural_filter(double(data_proj2(:,nt)), 1500, 0, 28, 32),  'r');
    
    % ==== plot cross spectru,
    j = J12(indf,nt);
    lt_subplot(nrow, 3, (3*nt)+3); hold on;
    plot(real(j), imag(j), 'ok');
    xlim([-30 30]);
    ylim([-30 30]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;


end



%% ======= 3) RECALCUALTE COHERENCE USING CROSS-TRIAL COHERENCE.









