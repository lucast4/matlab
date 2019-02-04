close all;
usecorrcoeff = 0; % for lfp xcorr/
% thingstodo = {'cohere'};
thintgstodo = {'lfpxcorr'};
% thingstodo = {'waveletcoh'};
[OUTSTRUCT, PARAMS, LFPXCORR_Base, LFPXCORR_WN, LFPXCORR_freqsall] =  ...
    lt_neural_LFP_RecalcCoh(OUTSTRUCT, SwitchCohStruct, LFPSTRUCT, SwitchStruct, ...
    PARAMS, usecorrcoeff, thingstodo);

%%

close all;
usecorrcoeff = 0; % for lfp xcorr/
% thingstodo = {'cohere'};
% thintgstodo = {'lfpxcorr'};
thingstodo = {'waveletcoh'};
[OUTSTRUCT, PARAMS, LFPXCORR_Base, LFPXCORR_WN, LFPXCORR_freqsall] =  ...
    lt_neural_LFP_RecalcCoh(OUTSTRUCT, SwitchCohStruct, LFPSTRUCT, SwitchStruct, ...
    PARAMS, usecorrcoeff, thingstodo);


%%

%% ====== EXTRACT BASELINE CORRELATIONS.
useonlybaseepoch = 0; % if 1, then just limited epoch. if 0, then entie baseline data
% default: 0;
corrtype = 'spearman'; % spearman. or pearson or Kendall
cohdiff_usedprime = 0; % if 0, then just mean diff. if 1, then dprime
nboot = 10; % to get bootstrap SE

OUTSTRUCT = lt_neural_Coher_CohScalExtract(OUTSTRUCT, useonlybaseepoch, ...
    corrtype, cohdiff_usedprime, nboot);



%%

%%

[OUTSTRUCT, OUTSTRUCT_CohMatOnly] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);


%% ################ COHERENCE SCALAR PLOTS [SYSTEMATIC]
corrtype = 'spearman'; % spearman. or pearson or Kendall
Nshuff = 500;

% GOES THRU ALL T, FF BINS
% REQUIRES OUTSTRUCT_CohMatOnly
tbins = PARAMS.tbins;
fbins = PARAMS.ffbins;

FFcorrCoh = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_pctileVsShuff = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_zscoreVsShuff = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_shuffCI = nan(length(tbins), length(fbins), 2, length(OUTSTRUCT.bnum));

savedir = '/bluejay5/lucas/analyses/neural/COHERENCE/SCALAR';
savedir = [savedir '/' PARAMS.savemarker];
if ~exist(savedir)
    mkdir(savedir);
end


for i=1:length(OUTSTRUCT.bnum)
    
    inds_base = OUTSTRUCT.indsbase{i};
    
    cohmatall = OUTSTRUCT_CohMatOnly{i};
    ff = OUTSTRUCT.ffvals{i};
    
    % ================ GO THRU EACH T, FF BIN
    for j=1:length(tbins)
        for jj=1:length(fbins)
            
            cohthis = squeeze(cohmatall(j, jj, inds_base));
            ffthis = ff(inds_base);
            
            % ---------------- GET CORR
            rho_dat = corr(ffthis', cohthis, 'type', corrtype);
            
            % ---------------- GET SHUFFLE DISTRIBUTION OF CORR
            rho_shuff_all = [];
            for nn=1:Nshuff
                disp(['case ' num2str(i) ', bin1=' num2str(j) ', bin2=' num2str(jj) ', shuff' num2str(nn)]);
                
                % --- shuffle trials
                ffthis_shuff = ffthis(randperm(length(ffthis)));
                rho_shuff = corr(ffthis_shuff', cohthis, 'type', corrtype);
                rho_shuff_all = [rho_shuff_all; rho_shuff];
                
            end
            
            % ================= OUTPUT STATS
            FFcorrCoh(j, jj, i) = rho_dat;
            
            p = (sum(abs(rho_shuff_all)>=abs(rho_dat))+1)./(length(rho_shuff_all)+1);
            FFcorrCoh_pctileVsShuff(j, jj, i) = p;
            
            shuffmean = mean(rho_shuff_all);
            shuffstd = std(rho_shuff_all);
            rho_z = (rho_dat - shuffmean)/shuffstd;
            FFcorrCoh_zscoreVsShuff(j, jj, i) = rho_z;
            FFcorrCoh_shuffCI(j, jj, :, i) = prctile(rho_shuff_all, [2.75 97.5]);
            
        end
    end
    
save([savedir '/FFcorrCoh'], 'FFcorrCoh');
save([savedir '/FFcorrCoh_pctileVsShuff'], 'FFcorrCoh_pctileVsShuff');
save([savedir '/FFcorrCoh_zscoreVsShuff'], 'FFcorrCoh_zscoreVsShuff');
save([savedir '/FFcorrCoh_shuffCI'], 'FFcorrCoh_shuffCI');

end


% clear OUTSTRUCT_CohMatOnly;