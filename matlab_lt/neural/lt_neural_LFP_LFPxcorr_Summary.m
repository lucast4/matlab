function lt_neural_LFP_LFPxcorr_Summary(OUTSTRUCT, PARAMS, SwitchStruct, ...
    birdstoplot, expttoplot, swtoplot, freqtoplot)
%% lt 12/15/18 - symmary of lfp xcorr analyses.
% CAn do any of the matrix data (e.g. S, coh...)

% sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs

assert(all(strcmp(OUTSTRUCT.bregionpair, 'LMAN-RA')), 'assumes chan1 is LAMN, chang 2 is RA...');
% tbins = PARAMS.tbins;
% ffbins = PARAMS.ffbins;

% clim = [-0.15 0.15];
plotEachTrial = 0; % for the timecourses, if 0 then only plots means, sem.

%% onluy plot certain birds?

OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT); % so followniog code works...

if ~isempty(birdstoplot)
    indstokeep = ismember(OUTSTRUCT.bnum, birdstoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

if ~isempty(expttoplot)
    indstokeep = ismember(OUTSTRUCT.enum, expttoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

if ~isempty(swtoplot)
    indstokeep = ismember(OUTSTRUCT.switch, swtoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end


%% ====== whic freq bin to plot?
indfreq = find(OUTSTRUCT.LFPXCORR_freqsall{1}==freqtoplot); assert(length(indfreq)==1);

OUTSTRUCT.LFPXCORR_FracChange = OUTSTRUCT.LFPXCORR_FracChange(:,indfreq);

%% ======= field type to plot
DATSTRUCT = struct;

% for coherence [diff]

% ============
%     fieldtoget = 'CohMean_WNminusBase';
%     [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
%         lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

fieldtoget = 'LFPXCORR_FracChange';
[~, ~, ~, ~, allbnum, allenum, allswnum, allDat] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

% fieldtoget = 'CohMean_Base';
% [~, ~, ~, ~, allbnum, allenum, allswnum, allDat1] = ...
%     lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
% 
% ---- subtract WN from base
% cohdiff = nan(size(allDat1));
% for i=1:size(cohdiff,3)
%     for ii=1:size(cohdiff,4)
%         
%         cohdiff(:,:,i,ii) = allDat2(:,:,i, ii) - allDat1(:,:,i, ii);
%     end
% end

% save dat
DATSTRUCT.lfpxcorr_diff = allDat;


%% ======================== PLOT targ, same, diff, for all switch

lt_figure; hold on;
ncases = size(allDat,4);
yall = [];
for j=1:ncases
   datthis = squeeze(allDat(:,:, :, j));
   
   x = 1:3;
   
   plot(x, datthis, '-ok');
   
   bname = SwitchStruct.bird(allbnum(j)).birdname;
   ename = SwitchStruct.bird(allbnum(j)).exptnum(allenum(j)).exptname;
   lt_plot_text(3.2, datthis(end), [bname '-' ename(end-8:end) '-sw' num2str(allswnum(j))], 'r');
   
   yall = [yall; datthis'];
end
% --- means
ymean = nanmean(yall,1);
ysem = lt_sem(yall);
lt_plot(x+0.2, ymean, {'Errors', ysem});

xlabel('TARG - SAME - DIFF');
ylabel('change in max(lfp-xcorr) [wn minus base]');
xlim([0 5]);
lt_plot_zeroline;



lt_figure; hold on;
ncases = size(allDat,4);
yall = [];
for j=1:ncases
   datthis = squeeze(allDat(:,:, [1 3], j));
   
   x = 1:2;
   
   plot(x, datthis, '-ok');
   
   bname = SwitchStruct.bird(allbnum(j)).birdname;
   ename = SwitchStruct.bird(allbnum(j)).exptnum(allenum(j)).exptname;
   lt_plot_text(2.2, datthis(end), [bname '-' ename(end-8:end) '-sw' num2str(allswnum(j))], 'r');
   
   yall = [yall; datthis'];
end
% --- means
ymean = nanmean(yall,1);
ysem = lt_sem(yall);
lt_plot(x+0.2, ymean, {'Errors', ysem});

xlabel('TARG - DIFF');
xlim([0 4]);
lt_plot_zeroline;

