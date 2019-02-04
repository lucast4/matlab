function lt_neural_Coher_COHSCALAR_Overview_sub(birdthis, OUTSTRUCT, FFcorrCoh, ...
    PARAMS, tlims, flims)

%% lt 1/8/19 - hierarhccail clustering based on correlation ofcorrleation (coh vs. ff) matrices.


%% get just this bird

indsgood = OUTSTRUCT.bnum==birdthis;
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indsgood, 1);
FFcorrCoh = FFcorrCoh(:,:, indsgood);

%% --------- RUN
inds_t = PARAMS.tbins>tlims(1) & PARAMS.tbins<tlims(2);
inds_ff = PARAMS.ffbins>flims(1) & PARAMS.ffbins<flims(2);

ncase = length(OUTSTRUCT.bnum);
distmat = nan(ncase, ncase);

XLAB = {};
for i=1:ncase
    disp(i);
    for ii=1:ncase
        
        cor1 = FFcorrCoh(inds_t, inds_ff, i);
        cor2 = FFcorrCoh(inds_t, inds_ff, ii);
        
        rho = corr(cor1(:), cor2(:));
        
        distmat(i, ii) = rho;
        
        
    end
    % -- GET LABEL
    % {motif - expt - sw}
    XLAB = [XLAB; ...
        ['mot' num2str(OUTSTRUCT.motifID_unique(i)) '-e' num2str(OUTSTRUCT.enum(i)) '-s' num2str(OUTSTRUCT.switch(i))]];
end



% ===================== HEIRARCHICAL CLUSTERING ON DISTANCE VECTORS
distmetric = 'euclidean';
% distmetric = 'correlation';
% locthis = 'RA';
% inds = find(strcmp(AllNeurLocation, locthis));

Z = linkage(distmat, 'single', distmetric);
Z = linkage(distmat, 'average', distmetric);
% lt_figure; hold on;
% dendrogram(Z, 0);
dendrogram(Z, 0, 'Labels', XLAB);
% dendrogram(Z, 'Labels', XLAB);
rotateXLabels(gca, 90);
% tmp = get(gca,'XTickLabel');


% ---- calcualte cophenetic correlation (i..e how well does clustering
% represent original distances?)

Y = pdist(distmat, distmetric);
c = cophenet(Z, Y);

%% =============== BREAK OUT ALL PAIRS IN DIFFERENT WAYS

% ==== 1) GET ALL PAIRWISE DISTANCES
distvec = pdist(distmat)';

% ==== SAME MOTIF VS. DIFF MOTIF


% ==== SAME EXPT VS. DIFF EXPT

