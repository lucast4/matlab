 %% ========== DATABASE OF BOS EXPERIMENTS
% --- ONE expt for each BOS presentation (at a specific depth)
clear all; close all;
ind = 0;


% % ========= NOTE, THIS IS A SMALL BATCH JUST FOR TESTING.
% ind = ind+1;
% SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
% SummaryBOS.expt(ind).batchname = 'Batch2353to2357';
% SummaryBOS.expt(ind).channels = [9 14 14];
% SummaryBOS.expt(ind).clusters = [1 1 2];
% SummaryBOS.expt(ind).bregions= {'test', 'test', 'test'};

% ================================== 11/9 - Night 0
% ind = ind+1; OLD VERSION. NOW EXCLUDING CHAN 8 SINCE IS NOT SONG MOD
% (DURING SINGING), SO NOT IN LMAN.
% SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
% SummaryBOS.expt(ind).batchname = 'Batch2327to2349good';
% SummaryBOS.expt(ind).channels = [8 9 9 14 21]; % note: 21 maybe coule be 2 clusteres?
% SummaryBOS.expt(ind).clusters = [1 1 2 1 1];
% SummaryBOS.expt(ind).isSU = [0 0 0 0 0];
% SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'LMAN', 'RA'};
ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
SummaryBOS.expt(ind).batchname = 'Batch2327to2349good';
SummaryBOS.expt(ind).channels = [9 9 14 21]; % note: 21 maybe coule be 2 clusteres?
SummaryBOS.expt(ind).clusters = [1 2 1 1];
SummaryBOS.expt(ind).isSU = [0 0 0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'RA'};

ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
SummaryBOS.expt(ind).batchname = 'Batch0004to0032';
SummaryBOS.expt(ind).channels = [9 14 14 17 21 21]; %
SummaryBOS.expt(ind).clusters = [1 1 2 1 1 2];
SummaryBOS.expt(ind).isSU = [0 0 1 0 0 1];
SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'RA', 'RA', 'RA'};


ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
SummaryBOS.expt(ind).batchname = 'Batch2353to2356';
SummaryBOS.expt(ind).channels = [9 14 14 21]; %
SummaryBOS.expt(ind).clusters = [1 1 2 1];
SummaryBOS.expt(ind).isSU = [0 0 0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'RA'};


% ================================== 11/10 - Night 1
ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night';
SummaryBOS.expt(ind).batchname = 'Batch0223to0251';
SummaryBOS.expt(ind).channels = [14 21]; %
SummaryBOS.expt(ind).clusters = [1 1];
SummaryBOS.expt(ind).isSU = [0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA'};


ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night';
SummaryBOS.expt(ind).batchname = 'Batch0158to0213';
SummaryBOS.expt(ind).channels = [14 21]; %
SummaryBOS.expt(ind).clusters = [1 1];
SummaryBOS.expt(ind).isSU = [0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA'};


% ================================== 11/11 - Night 2
ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111117Night';
SummaryBOS.expt(ind).batchname = 'Batch0202to0226';
SummaryBOS.expt(ind).channels = [14 21 21]; %
SummaryBOS.expt(ind).clusters = [1 1 2];
SummaryBOS.expt(ind).isSU = [0 0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA', 'RA'};

ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111117Night';
SummaryBOS.expt(ind).batchname = 'Batch0232to0241';
SummaryBOS.expt(ind).channels = [14 21 21]; %
SummaryBOS.expt(ind).clusters = [1 1 2];
SummaryBOS.expt(ind).isSU = [0 0 1];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA', 'RA'};


% ####################### AUTO EXTRACTION OF BIRD NAME
for i=1:length(SummaryBOS.expt)
    tmp1 = strfind(SummaryBOS.expt(i).dirname, 'birds/');
    tmp2 = strfind(SummaryBOS.expt(i).dirname, '/NEURAL');
    bname = SummaryBOS.expt(i).dirname(tmp1+6:tmp2-1);
    SummaryBOS.expt(i).birdname = bname;
end

% ======== sanity check (entered one thing for each unit..)
for i=1:length(SummaryBOS.expt)
    assert(length(unique([length(SummaryBOS.expt(i).channels); ...
        length(SummaryBOS.expt(i).clusters); ...
        length(SummaryBOS.expt(i).isSU); ...
        length(SummaryBOS.expt(i).bregions); ...
        ]))==1, 'did not enter correcly, all shoudl have same number units');
end

%% ##################### [MAIN EXTRACTION]
close all;

expttoget = 1;
assert(length(SummaryBOS.expt) ==1, 'currently only sure this works well if only have one expt loaded..');
lt_neural_BOS_Extraction(SummaryBOS, expttoget)


% ============= NOTE: following this, cheeck using wave_clus.

%% ============ 4) PLOT EACH SONG FILE OVERLAID WITH EXTRACTED SPIKES
i = 1;
fs = 30000;
batchname = SummaryBOS.expt(i).batchname;
chanstoget = 14;

% ======================= GO THRU EACH FILENAME
filenames = textread(batchname, '%s');
figcount=1;
subplotrows=3;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for j=1:length(filenames)
    fname = filenames{j};
    
    % ============= load and extract each channel
    % ---- 1) raw neural
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fname);
    
    [a, b] = fileparts(fname);
    
    % === go thru each channle
    for jj=1:length(chanstoget)
        chan = chanstoget(jj);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([fname '-chan' num2str(chan)]);
        
        neurdat = load([batchname '-Chan' num2str(chan) '/' b '.ch' num2str(chan) '.mat']);
        neurdatold = amplifier_data([amplifier_channels.chip_channel] == chan,:);
        assert(all(neurdat.data==single(neurdatold)));
        
        % --- filter neural dat
        neurdat = lt_neural_filter(double(neurdat.data), fs, 1);
        
        % ------ EXTRACT SPIKES
        %         spksold = load([batchname '-Chan' num2str(chan) '/' b '.ch' num2str(chan) '_spikes.mat']);
        spksnew = load([batchname '-Chan' num2str(chan) '/allsongs_spikes.mat']);
        timesdat = load([batchname '-Chan' num2str(chan) '/times_allsongs.mat']);
        %         spksnew.filenuminbatch_all==j;
        assert(size(spksnew.index,2) == size(timesdat.cluster_class,1));
        
        clustclassthis = timesdat.cluster_class(spksnew.filenuminbatch_all==j, :);
        
        % ================================ PLOT STUFF
        t = 1:length(neurdat);
        t = (t./fs)*1000;
        plot(t, neurdat, 'k');
        numclust = max(clustclassthis(:,1));
        plotcols = lt_make_plot_colors(numclust, 0, 0);
        for jjj=1:numclust % go thru all clusters (ignore 0, which is noise)
            spkthis = clustclassthis(clustclassthis(:,1)==jjj,2);
            lt_plot(spkthis, jjj, 'o', 'Color', plotcols{jjj});
        end
    end
end



%% ##############################################################
%% ############################ [ANALYSES]
%% ========= PLOT OVERVIEW OF DATA
i = 1;

ls(SummaryBOS.expt(i).dirname, '-l');
disp('LIST OF ALL BOS SONGS');


%% ======================= [EXTRACTION]
% ======== FIRST RUN THIS TO EXTRACT ALL PREVIOUSLY PROCESSED DATA.

lt_neural_BOS_Extraction_Script;

%% ##################### [PLOT] SUMMARY FOR RESPONSE [BOSTYPES] - song level
% EACH UNIT PLOT MEAN RESPONSE TO EACH BOS TYPE
close all;
i=1;
birdthis = SummaryBOS.expt(i).birdname;

lt_neural_BOS_SummaryPlot;


%% ###################### [MOTIF-LEVEL PLOTS]

lt_neural_BOS_MotifScript;



%% ########################## [COHERENCE]
for i=1:length(SummaryBOS.expt)
    
    % ============================== PARAMS
    lt_switch_chronux(1);
    movingwin = [0.1 0.01];
    params = struct;
    params.fpass = [1/movingwin(1) 200];
    w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
    % w = 20; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
    tw = movingwin(1)*w;
    params.tapers = [tw 2*tw-1];
    params.Fs = 1500; % hard coded fs for LFP;
    
    
    % =============================
    LFPall = SummaryBOS.expt(i).DAT_bysongrend.LFPall;
    LFPchanlist = SummaryBOS.expt(i).DAT_bysongrend.LFPall_chanlist;
    LFP_t = SummaryBOS.expt(i).DAT_bysongrend.LFPall_t;
    
    nsongs = size(LFPall,1);
    nchans = unique(cellfun(@length, SummaryBOS.expt(i).DAT_bysongrend.LFPall_chanlist));
    
    CohCell = cell(nchoosek(nchans, 2), nsongs); % pairs x songs
    Coh_t =  cell(1, nsongs); % pairs x songs
    Coh_f = cell(1, nsongs); % pairs x songs
    Coh_chanpair = cell(nchoosek(nchans, 2), nsongs); % pairs x songs
    Coh_bregionpair = cell(nchoosek(nchans, 2), nsongs); % pairs x songs\
    
    chpair = [];
    for ss=1:nsongs
        count = 1; % over chan p[aiors.
        disp(['song ' num2str(ss)]);
        for cc=1:nchans
            for ccc=cc+1:nchans
                
                lfp1 = LFPall{ss, cc};
                lfp2 = LFPall{ss, ccc};
                
                tthis = LFP_t{ss};
                
                chan1 = LFPchanlist{ss}(cc);
                chan2 = LFPchanlist{ss}(ccc);
                
                % --- bregions
                bregion1 = unique(SummaryBOS.expt(i).bregions(SummaryBOS.expt(i).channels == chan1));
                bregion2 = unique(SummaryBOS.expt(i).bregions(SummaryBOS.expt(i).channels == chan2));
                bregion1 = bregion1{1};
                bregion2 = bregion2{1};
                
                if strcmp(bregion1, 'RA')
                    assert(~strcmp(bregion2, 'LMAN'), 'if LMAN and RA, then must be LMAN first [assumed later on]');
                end
                
                
                % ===== get coherogram over entire song
                [C,~,~,~,~,t,f] = cohgramc(lfp1, lfp2, movingwin, params);
                %             phi = single(phi);
                %             S12 = single(S12);
                %             S1 = single(S1);
                %             S2 = single(S2);
                
                % --- convert t to relative to dat onset
                t = t+PARAMS.flanktotake(1);
                
                % ====== SAVE COHEROGRAM
                CohCell{count, ss} = single(C);
                Coh_chanpair{count, ss} = [chan1 chan2];
                Coh_bregionpair{count, ss} = [bregion1 '-' bregion2];
                
                if count==1
                    Coh_t{ss} = [single(t(1)) single(t(end))];
                    Coh_f{ss} = single(f);
                end
                
                count = count+1;
                
                
                
                % =========================================
                if (0) % to plot
                    figure; hold on; lt_neural_Coher_Plot(C, t, f, 1, '', [], 0, 0);
                    lt_plot_colormap('pval');
                    
                    %                line(PARAMS.flanktotake
                    
                end
            end
            
        end
    end
    
    lt_switch_chronux(0);
    
    % ==== SAVE OUTPUT
    SummaryBOS.expt(i).DAT_bysongrend.CohCell = CohCell';
    SummaryBOS.expt(i).DAT_bysongrend.Coh_chanpair = Coh_chanpair';
    SummaryBOS.expt(i).DAT_bysongrend.Coh_bregionpair = Coh_bregionpair';
    SummaryBOS.expt(i).DAT_bysongrend.Coh_t = Coh_t';
    SummaryBOS.expt(i).DAT_bysongrend.Coh_f = Coh_f';
end

%% ============ [COPHERENCE] - PLOT ALL
close all;
numexpts = length(SummaryBOS.expt);
bostypelist = {'fwd', 'rev'}; % if is length 2, will also get difference
bregiontoget = 'LMAN-RA';
clim = [0.2 0.9];
windtoplot_base = [-2.5 -1]; % rel onset
windtoplot_BOS = [1 10]; % rel onset
fwindtoplot = [20 32];
onlyplotFirstSong = 1; % i.e. for each expt, only the first trial. (e.g. to avoid habituation)


figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


Y = []; % base, bos
for i=1:numexpts
    
    % -- CHANLLE PAIRS
    tmp = cellfun(@(x)strcmp(x, bregiontoget), SummaryBOS.expt(i).DAT_bysongrend.Coh_bregionpair);
    for j=1:size(tmp,2)
        assert(length(unique(tmp(:,j)))==1, 'why some trials diff brain regions?')
    end
    chanpairstoget_all = find(tmp(1,:));
    
    for cc=1:length(chanpairstoget_all)
        chanthis = chanpairstoget_all(cc);
        
        for ii=1:length(bostypelist)
            
            bostoget = bostypelist{ii};
            birdnum = strcmp(PARAMS.BOSbirdname, SummaryBOS.expt(i).birdname);
            bostoget_ind = find(strcmp(PARAMS.BOSnames{birdnum}, bostoget));
            
            
            % ===== what trials and chanpairs to get
            % -- TRIALS
            inds_trials = SummaryBOS.expt(i).DAT_bysongrend.BOStype==bostoget_ind;
%             % -- CHANLLE PAIRS
%             tmp = cellfun(@(x)strcmp(x, bregiontoget), SummaryBOS.expt(i).DAT_bysongrend.Coh_bregionpair);
%             for j=1:size(tmp,2)
%                 assert(length(unique(tmp(:,j)))==1, 'why some trials diff brain regions?')
%             end
%             chanpairstoget = find(tmp(1,:));
            if onlyplotFirstSong==1
                inds_trials = intersect(find(SummaryBOS.expt(i).DAT_bysongrend.BOStype==bostoget_ind), ...
                    1);
            end
              
            if isempty(inds_trials)
                continue
            end
               
            
            % ========== COLLECT DAT
            datall = SummaryBOS.expt(i).DAT_bysongrend.CohCell(inds_trials, chanthis);
            datall = datall(:);
            datall = lt_neural_Coher_Cell2Mat(datall);
            
            t = linspace(SummaryBOS.expt(i).DAT_bysongrend.Coh_t{1}(1), ...
                SummaryBOS.expt(i).DAT_bysongrend.Coh_t{1}(2), size(datall,1));
            f = SummaryBOS.expt(i).DAT_bysongrend.Coh_f{1};
            
            % ============= PLOT 1 - COHEROGRAM
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['expt' num2str(i) '[' bostoget '](' bregiontoget ')' ]);
            lt_neural_Coher_Plot(datall, t, f, 1, '', clim);
            lt_plot_colormap('pval');
            lt_plot_annotation(1, ['n=' num2str(size(datall,3))], 'm')
            
            % ============= PLOT 2 -
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['expt' num2str(i) '[' bostoget '](' bregiontoget ')' ]);
            lt_neural_Coher_Plot(datall, t, f, 2, '-', [0 1], 0, 0, fwindtoplot);
            %         lt_plot_colormap('pval');
           
            % ======================== PLOT MEAN COH SPECTRUM
            ind_base = t>windtoplot_base(1) & t<windtoplot_base(2);
            ind_BOS = t>windtoplot_BOS(1) & t<windtoplot_BOS(2);
            
            cohspec_base = mean(mean(datall(ind_base, :, :),1),3);
            cohspec_base_sem = lt_sem(squeeze(mean(datall(ind_base, :, :),1))');
            cohspec_BOS = mean(mean(datall(ind_BOS, :, :),1),3);
            cohspec_BOS_sem = lt_sem(squeeze(mean(datall(ind_BOS, :, :),1))');
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            ylabel('coh');
            title('k=base; rd=BOS');
            xlabel('f(hz)');
            if length(cohspec_base_sem)>1
            shadedErrorBar(f, cohspec_base, cohspec_base_sem, {'Color', 'k'},1);
            shadedErrorBar(f, cohspec_BOS, cohspec_BOS_sem, {'Color', 'r'},1);
            end
        axis tight;

        end
        
        % ======== PLOT DIFFERENCE
        % ----- 1) FWD
        bostoget = 'fwd';
        
        birdnum = strcmp(PARAMS.BOSbirdname, SummaryBOS.expt(i).birdname);
        bostoget_ind = find(strcmp(PARAMS.BOSnames{birdnum}, bostoget));
        
        
        % ===== what trials and chanpairs to get
        % -- TRIALS
        inds_trials = SummaryBOS.expt(i).DAT_bysongrend.BOStype==bostoget_ind;
        %     % -- CHANLLE PAIRS
        %     tmp = cellfun(@(x)strcmp(x, bregiontoget), SummaryBOS.expt(i).DAT_bysongrend.Coh_bregionpair);
        %     for j=1:size(tmp,2)
        %         assert(length(unique(tmp(:,j)))==1, 'why some trials diff brain regions?')
        %     end
        %     chanpairstoget = find(tmp(1,:));
        
        % ========== COLLECT DAT
        datall = SummaryBOS.expt(i).DAT_bysongrend.CohCell(inds_trials, chanthis);
        datall = datall(:);
        datall = lt_neural_Coher_Cell2Mat(datall);
        
        % --
        datall_FWD = datall;
        
        
        
        % ----- 1) REV
        bostoget = 'rev';
        
        birdnum = strcmp(PARAMS.BOSbirdname, SummaryBOS.expt(i).birdname);
        bostoget_ind = find(strcmp(PARAMS.BOSnames{birdnum}, bostoget));
        
        
        % ===== what trials and chanpairs to get
        % -- TRIALS
        inds_trials = SummaryBOS.expt(i).DAT_bysongrend.BOStype==bostoget_ind;
        %     % -- CHANLLE PAIRS
        %     tmp = cellfun(@(x)strcmp(x, bregiontoget), SummaryBOS.expt(i).DAT_bysongrend.Coh_bregionpair);
        %     for j=1:size(tmp,2)
        %         assert(length(unique(tmp(:,j)))==1, 'why some trials diff brain regions?')
        %     end
        %     chanpairstoget = find(tmp(1,:));
        
        % ========== COLLECT DAT
        datall = SummaryBOS.expt(i).DAT_bysongrend.CohCell(inds_trials, chanthis);
        datall = datall(:);
        datall = lt_neural_Coher_Cell2Mat(datall);
        
        % --
        datall_REV = datall;
        
        
        % ======================================== DIFFERENCES
        datall = mean(datall_FWD,3) - mean(datall_REV,3);
        
        % ============= PLOT 1 - COHEROGRAM
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['expt' num2str(i) '[fwd-rev](' bregiontoget ')' ]);
        lt_neural_Coher_Plot(datall, t, f, 1, '', [-0.2 0.2]);
        lt_plot_colormap('centered');
        % -- overlay windows
        line([windtoplot_base(1) windtoplot_base(1)], ylim, 'Color', 'y');
        line([windtoplot_base(2) windtoplot_base(2)], ylim, 'Color', 'y');
        line([windtoplot_BOS(1) windtoplot_BOS(1)], ylim, 'Color', 'w');
        line([windtoplot_BOS(2) windtoplot_BOS(2)], ylim, 'Color', 'w');
        
        % ============= PLOT 2 -
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %     title(['expt' num2str(i) '[' bostoget '](' bregiontoget ')' ]);
        lt_neural_Coher_Plot(datall, t, f, 2, '-', [-0.4 0.4], 1, 0, fwindtoplot);
        %         lt_plot_colormap('pval');
        lt_plot_zeroline;
        

        %% ============= collapsing across time [and freq sometimes]
       
        % ======================== PLOT MEAN COH SPECTRUM
        ind_base = t>windtoplot_base(1) & t<windtoplot_base(2);
        ind_BOS = t>windtoplot_BOS(1) & t<windtoplot_BOS(2);
        ind_f = f>fwindtoplot(1) & f<fwindtoplot(2);
        
        cohspec_base = mean(datall(ind_base, :),1);
        cohspec_base_sem = lt_sem(datall(ind_base, :));
        cohspec_BOS = mean(datall(ind_BOS, :),1);
        cohspec_BOS_sem = lt_sem(datall(ind_BOS, :));
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('coh (fwd - rev)');
        title('k=base; rd=BOS');
        xlabel('f(hz)');
        shadedErrorBar(f, cohspec_base, cohspec_base_sem, {'Color', 'k'},1);
        shadedErrorBar(f, cohspec_BOS, cohspec_BOS_sem, {'Color', 'r'},1);
        
        % ======================== COLLECT DIFFERENCES, SCALAR
        if size(datall_REV,3)<5 | size(datall_FWD,3)<5
            continue
        end
         y1 = mean(mean(datall(ind_base, :),1),2);
        y2 = mean(mean(datall(ind_BOS, :),1),2);
        
       
        Y = [Y; [y1 y2]];
    end
end

lt_figure; hold on;
plot([1 2], Y', '-ko');
title('each chan pair');
xlabel('baseline -- BOS')
ylabel('fwd minus rev');
xlim([0 3]);

%% ====================== [PLOT] RASTERS OVER TRIALS, ONE FOR EACH UNIT
% NOTE: quite slow since plots a lot of rasters.

close all;
i = 2;
lt_neural_BOS_Rasters;



%% +++++++++++++++++++++++++++++++++++++++++ OLD STUFF

%% ====================== [OLD] OLD ANALYSES - IGNORE GENRALLY

lt_neural_BOS_OLDSTUFF;

%% ######################################################################
%% ######################################################################
%%  TO MAKE BOS FILES (with digitual pulse for syl onset/offset

lt_neural_BOS_MakeBOS_Script;
