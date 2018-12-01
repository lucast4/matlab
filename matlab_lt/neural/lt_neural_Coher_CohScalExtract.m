function OUTSTRUCT = lt_neural_Coher_CohScalExtract(OUTSTRUCT, useonlybaseepoch, ...
    corrtype, cohdiff_usedprime, nboot)

% useonlybaseepoch = 0; % if 1, then just limited epoch. if 0, then entie baseline data
% % default: 0;
% corrtype = 'spearman'; % spearman. or pearson or Kendall
% cohdiff_usedprime = 0; % if 0, then just mean diff. if 1, then dprime
% nboot = 10; % to get bootstrap SE

%%
OUTSTRUCT.cohscal_ff_rho_rhoSE_base = [];
OUTSTRUCT.cohscal_WN_minus_base = [];

%%
% ============= 1) get baseline correlation between coh and FF
for i=1:length(OUTSTRUCT.bnum)
   
    disp(num2str(i));
    if useonlybaseepoch==1
        inds_base = OUTSTRUCT.indsbase_epoch{i};
    else
        inds_base = OUTSTRUCT.indsbase{i};
    end
    
    cohscal = OUTSTRUCT.cohscal{i}(inds_base);
    ff = OUTSTRUCT.ffvals{i}(inds_base);
    tvals = OUTSTRUCT.tvals{i}(inds_base);
    
    % ============ baseline corr
    Qfunc = @(x)corr(x(:,1), x(:,2), 'type', corrtype);
    X = [ff', cohscal'];
    rho_base = Qfunc(X);
    if ~isnan(rho_base)
        assert(corr(ff', cohscal', 'type', corrtype)==rho_base);
    end
    % --- get bootstrap 
    if isnan(rho_base)
        rho_SE = nan;
    else
        rho_shuff_all = [];
        for nn=1:nboot
            indshuff = randi(length(ff), length(ff), 1);
            x_shuff = X(indshuff, :);
            
            rho_shuff_all = [rho_shuff_all; Qfunc(x_shuff)];
        end
        rho_SE = std(rho_shuff_all);
    end
    
    
    % ============ WN minus base, change in coh
    cohscal = OUTSTRUCT.cohscal{i};
    indsWN = OUTSTRUCT.indsWN_epoch{i};
    indsbase = OUTSTRUCT.indsbase_epoch{i};
    if cohdiff_usedprime==0
    cohscal_WnMinusBase = mean(cohscal(indsWN) - mean(cohscal(indsbase)));
    else
        asdfasdf;
    end
    
    % ============== SAVE
    OUTSTRUCT.cohscal_ff_rho_rhoSE_base = [OUTSTRUCT.cohscal_ff_rho_rhoSE_base; ...
        [rho_base rho_SE]];
    OUTSTRUCT.cohscal_WN_minus_base = [OUTSTRUCT.cohscal_WN_minus_base; ...
        cohscal_WnMinusBase];
    
end