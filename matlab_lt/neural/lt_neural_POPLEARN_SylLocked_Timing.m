function lt_neural_POPLEARN_SylLocked_Coord(DATSTRUCT_SPK, DATSTRUCT_LFP, ...
    PARAMS, plotRawRand)
%% lt 2/25/19 - 

%% ===
% disp('NOTE: any cases with nan for correlation I replace by median acros tirlas (for that neuron pair x syllable)');

% =========== PARAMS FOR COHERENCE
% ntapers = []; % leave [] to set as default.
% movingwin = [0.1 0.01]; % leave exmpty for default. [applies for both welches and multiutaper]
% tw = [];
%
% fwind = [22 36]; % window to get scalar;

% ==== match movingwin to size of data
% wintmp = size(DATSTRUCT_LFP.LFP_dat{1},1)/1500;
% disp(['changin cohernece window from: ' num2str(movingwin) ', to: ' num2str([wintmp movingwin(2)])]);
% movingwin(1) = wintmp;


% =========== cross correlations
% corrtype = 'coeff'; % if this, then will plot the corr coeff
corrtype = 'unbiased'; % if this, then will z-score the final trial-averaged xcorr so that units are matched across neuron pairs
maxlag = 0.05; % seconds;

%% ============= for each experiment, get coordination in spiking
[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_SPK.bnum, DATSTRUCT_SPK.enum, ...
    DATSTRUCT_SPK.switch, DATSTRUCT_SPK.motifnum});

DATSTRUCT_POP = struct;
DATSTRUCT_POP.bnum = [];
DATSTRUCT_POP.enum = [];
DATSTRUCT_POP.switch = [];
DATSTRUCT_POP.motifnum = [];
DATSTRUCT_POP.bregion = {};

DATSTRUCT_POP.xcorr_lfp = {};
DATSTRUCT_POP.xcorr_frsmth = {};

for i=1:length(indsgrpU)
    disp(i);
    
    indsthis = find(indsgrp==indsgrpU(i));
    
    bnum = unique(DATSTRUCT_SPK.bnum(indsthis));
    enum = unique(DATSTRUCT_SPK.enum(indsthis));
    sw = unique(DATSTRUCT_SPK.switch(indsthis));
    mm = unique(DATSTRUCT_SPK.motifnum(indsthis));
    
    % ######################################### LMAN
    bregionthis = 'LMAN';
    
    indsthis_spk = find(indsgrp==indsgrpU(i) & ...
        strcmp(DATSTRUCT_SPK.spike_bregions, bregionthis));
    
    indsthis_lfp = DATSTRUCT_SPK.LFP_indthis(indsthis_spk);
    
    % --------------- COLLECT DATA
    indstrials = DATSTRUCT_SPK.inds_base{indsthis_spk(1)};
    frmat = DATSTRUCT_SPK.frmat(indsthis_spk);
    spkdat = DATSTRUCT_SPK.spike_dat(indsthis_spk);
    
    lfpdat = DATSTRUCT_LFP.LFP_dat(indsthis_lfp);
    
    assert(size(frmat{1},2) == size(lfpdat{1},2));
    
    % ---------- GET TRIALS OF INTEREST
    frmat = cellfun(@(x)x(:, indstrials), frmat, 'UniformOutput', 0);
%     spkdat = cellfun(@(x)x(indstrials), spkdat, 'UniformOutput', 0);
    lfpdat = cellfun(@(x)x(:, indstrials), lfpdat, 'UniformOutput', 0);
    
    % -- save for later
    lfpdat_LMAN = lfpdat;
    %     spkdat_LMAN = spkdat;
    frmat_LMAN = frmat;
    
    % ######################################### RA
    bregionthis = 'RA';
    
    indsthis_spk = find(indsgrp==indsgrpU(i) & ...
        strcmp(DATSTRUCT_SPK.spike_bregions, bregionthis));
    
    indsthis_lfp = DATSTRUCT_SPK.LFP_indthis(indsthis_spk);
    
    % --------------- COLLECT DATA
    indstrials = DATSTRUCT_SPK.inds_base{indsthis_spk(1)};
    frmat = DATSTRUCT_SPK.frmat(indsthis_spk);
    spkdat = DATSTRUCT_SPK.spike_dat(indsthis_spk);
    
    lfpdat = DATSTRUCT_LFP.LFP_dat(indsthis_lfp);
    
    assert(size(frmat{1},2) == size(lfpdat{1},2));
    
    % ---------- GET TRIALS OF INTEREST
    frmat = cellfun(@(x)x(:, indstrials), frmat, 'UniformOutput', 0);
%     spkdat = cellfun(@(x)x(indstrials), spkdat, 'UniformOutput', 0);
    lfpdat = cellfun(@(x)x(:, indstrials), lfpdat, 'UniformOutput', 0);
    
    % -- save for later
    lfpdat_RA = lfpdat;
    %     spkdat_RA = spkdat;
    frmat_RA = frmat;
    
    
    
    % ##################################### GET ALL CROSS CORRELATIONS
    % ============================== SPIKE VS. SPIKE
    t_fr = PARAMS.THIS.frmat_x;
    binsize = t_fr(2)-t_fr(1);
    maxlagthis = round(maxlag/binsize);
    ntrials = size(frmat_LMAN{1},2);
    
    nchan_LMAN = length(frmat_LMAN);
    nchan_RA = length(frmat_RA);
    
    % --- go thru all pairs of LMAN and RA channels, once for each trials
    ccorr_all = cell(1, nchan_LMAN*nchan_RA);
    cc=1;
    for nn=1:nchan_LMAN
        for nnn=1:nchan_RA
            
            xcorrall = [];
            for tt=1:ntrials
                
                frL = frmat_LMAN{nn}(:, tt);
                frR = frmat_RA{nnn}(:, tt);
                
                [y, lags1] = xcorr(frL, frR, maxlagthis, corrtype);
                
                %                             [y2, lags] = xcov(frL, frR, maxlagthis, corrtype);
                %                             [y2, lags] = xcorr(frL, frR, maxlagthis, 'unbiased');
                
                
                % ========= SAVE
                xcorrall = [xcorrall; y'];
                
                %                             if (0)
                %                                figure
                %
                %                                lt_subplot(3,2,1); hold on;
                %                                title('k=LMAN');
                %                                plot(t_fr, frL, 'k');
                %                                plot(t_fr, frR, '-r');
                %
                %                                lt_subplot(3,2,2); hold on;
                %                                xlabel('LM<AN-RA');
                %                                title('xcorr(r), xcov (b)');
                %                                plot(lags1, y, 'r');
                %                                plot(lags, y2, 'b');
                %
                %                                lt_subplot(3,2,3); hold on;
                %                                xlabel('LM<AN-RA');
                %                                title('xcorr(r), xcov (b)');
                %                 %                plot(lags1, y, 'r');
                %                                plot(lags, y2, 'b');
                %                             end
            end
            
            ccorr_all{cc} = xcorrall';
            cc = cc+1;
            PARAMS.THIS.xcorr_fr_lags = lags1*binsize;
        end
    end
    DATSTRUCT_POP.xcorr_frsmth = [DATSTRUCT_POP.xcorr_frsmth; {ccorr_all}];
    
    
    % =============================== LFP VS. LFP
    t_LFP = PARAMS.THIS.lfpx;
    binsize = t_LFP(2)-t_LFP(1);
    maxlagthis = round(maxlag/binsize);
    ntrials = size(lfpdat_LMAN{1},2);
    nchan_LMAN = length(lfpdat_LMAN);
    nchan_RA = length(lfpdat_RA);
    
    % --- go thru all pairs of LMAN and RA channels, once for each trials
    ccorr_all = cell(1, nchan_LMAN*nchan_RA);
    cc=1;
    for nn=1:nchan_LMAN
        for nnn=1:nchan_RA
            
            xcorrall = [];
            for tt=1:ntrials
                
                lfpL = lfpdat_LMAN{nn}(:, tt);
                lfpR = lfpdat_RA{nnn}(:, tt);
                
                [y, lags1] = xcorr(lfpL, lfpR, maxlagthis, corrtype);
                
                %             [y2, lags] = xcov(lfpL, lfpR, [], corrtype);
                %             [y2, lags] = xcorr(lfpL, lfpR, [], 'unbiased');
                
                
                % ========= SAVE
                xcorrall = [xcorrall; y'];
                
                %             if (0)
                %                figure
                %
                %                lt_subplot(3,2,1); hold on;
                %                title('k=LMAN');
                %                plot(t_LFP, lfpL, 'k');
                %                plot(t_LFP, lfpR, '-r');
                %
                %                lt_subplot(3,2,2); hold on;
                %                xlabel('LM<AN-RA');
                %                title('xcorr(r), xcov (b)');
                %                plot(lags1, y, 'r');
                %                plot(lags, y2, 'b');
                %
                %                lt_subplot(3,2,3); hold on;
                %                xlabel('LM<AN-RA');
                %                title('xcorr(r), xcov (b)');
                % %                plot(lags1, y, 'r');
                %                plot(lags, y2, 'b');
                %             end
            end
            
            ccorr_all{cc} = xcorrall';
            cc = cc+1;
            PARAMS.THIS.xcorr_lfp_lags = lags1*binsize;
        end
    end
    DATSTRUCT_POP.xcorr_lfp = [DATSTRUCT_POP.xcorr_lfp; {ccorr_all}];
    
    
    % ===================== ADD OTHER THINGS    
    DATSTRUCT_POP.bnum = [DATSTRUCT_POP.bnum; bnum];
    DATSTRUCT_POP.enum = [DATSTRUCT_POP.enum; enum];
    DATSTRUCT_POP.switch = [DATSTRUCT_POP.switch; sw];
    DATSTRUCT_POP.motifnum = [DATSTRUCT_POP.motifnum; mm];
    DATSTRUCT_POP.bregion = [DATSTRUCT_POP.bregion; 'LMAN-RA'];
end


%% ===== [PLOT]
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% -- 
if strcmp(corrtype, 'unbiased')
    takezscoreofmean =1; % if 1, then zscores the xcorr trace [acuitl;al;y subtracts mean only!]
else
    takezscoreofmean =0; % if 1, then zscores the xcorr trace
end

% ==================== 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('trial-by-trial xcorr (coeff');

y = DATSTRUCT_POP.xcorr_lfp;
x = PARAMS.THIS.xcorr_lfp_lags;
% --- for each case get trial-average
for i=1:length(y)
    y{i} = cellfun(@(x)mean(x,2), y{i}, 'UniformOutput', 0);
end
y = cellfun(@(x)cell2mat(x), y, 'UniformOutput', 0); % concatenate pairs
y = cellfun(@(x)x', y, 'UniformOutput', 0);
y = cell2mat(y)';
% - take zscore for each case
if takezscoreofmean==1
    y = (y - mean(y,1))./std(y);
%     y = (y - mean(y,1));
end
% -- sort by max
[~, indmax] = max(y);
[~, indsort] = sort(indmax);
y = y(:, indsort);
clim = [min(y(:)) max(y(:))];
f = 1:size(y,2);
lt_neural_Coher_Plot(y, x, f, 1, '', clim);


% ==================== 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('trial-by-trial xcorr (coeff');
xlabel('LMAN - RA');
ylabel('LFP-LFP');

y = DATSTRUCT_POP.xcorr_lfp;
x = PARAMS.THIS.xcorr_lfp_lags;
% --- for each case get trial-average
for i=1:length(y)
    y{i} = cellfun(@(x)mean(x,2), y{i}, 'UniformOutput', 0);
end
y = cellfun(@(x)cell2mat(x), y, 'UniformOutput', 0); % concatenate pairs
y = cellfun(@(x)x', y, 'UniformOutput', 0);
y = cell2mat(y)';
% - take zscore for each case
if takezscoreofmean==1
    y = (y - mean(y,1))./std(y);
%     y = (y - mean(y,1));
end

ymean = mean(y,2);
ysem = lt_sem(y');
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);


% ==================== 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('trial-by-trial xcorr (coeff');

y = DATSTRUCT_POP.xcorr_frsmth;
x = PARAMS.THIS.xcorr_fr_lags;
% --- for each case get trial-average
for i=1:length(y)
    y{i} = cellfun(@(x)mean(x,2), y{i}, 'UniformOutput', 0);
end
y = cellfun(@(x)cell2mat(x), y, 'UniformOutput', 0); % concatenate pairs
y = cellfun(@(x)x', y, 'UniformOutput', 0);
y = cell2mat(y)';
% - take zscore for each case
if takezscoreofmean==1
    y = (y - mean(y,1))./std(y);
%     y = (y - mean(y,1));
end
% -- sort by max
[~, indmax] = max(y);
[~, indsort] = sort(indmax);
y = y(:, indsort);
clim = [min(y(:)) max(y(:))];
f = 1:size(y,2);
lt_neural_Coher_Plot(y, x, f, 1, '', clim);


% ==================== 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('trial-by-trial xcorr (coeff');
xlabel('LMAN - RA');
ylabel('SPK-SPK');

y = DATSTRUCT_POP.xcorr_frsmth;
x = PARAMS.THIS.xcorr_fr_lags;
% --- for each case get trial-average
for i=1:length(y)
    y{i} = cellfun(@(x)mean(x,2), y{i}, 'UniformOutput', 0);
end
y = cellfun(@(x)cell2mat(x), y, 'UniformOutput', 0); % concatenate pairs
y = cellfun(@(x)x', y, 'UniformOutput', 0);
y = cell2mat(y)';
% - take zscore for each case
if takezscoreofmean==1
    y = (y - mean(y,1))./std(y);
%     y = (y - mean(y,1));
end
ymean = nanmean(y,2);
ysem = lt_sem(y');
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

