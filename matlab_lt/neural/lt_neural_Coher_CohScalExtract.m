function OUTSTRUCT = lt_neural_Coher_CohScalExtract(OUTSTRUCT, useonlybaseepoch, ...
    corrtype, cohdiff_usedprime, nboot)

% useonlybaseepoch = 0; % if 1, then just limited epoch. if 0, then entie baseline data
% % default: 0;
% corrtype = 'spearman'; % spearman. or pearson or Kendall
% cohdiff_usedprime = 0; % if 0, then just mean diff. if 1, then dprime
% nboot = 10; % to get bootstrap SE

%%
OUTSTRUCT.cohscal_ff_rho_rhoSE_base = [];
OUTSTRUCT.cohscal_ff_rho_rhoSE_WN= [];
OUTSTRUCT.cohscal_ff_rho_rhoSE_All= [];

OUTSTRUCT.cohscal_WN_minus_base = [];

%%
% ============= 1) get baseline correlation between coh and FF
for i=1:length(OUTSTRUCT.bnum)
    
    disp(num2str(i));
    %     tvals = OUTSTRUCT.tvals{i}(inds_base);
    
    %% ============ baseline corr
    if useonlybaseepoch==1
        indsthis = OUTSTRUCT.indsbase_epoch{i};
    else
        indsthis = OUTSTRUCT.indsbase{i};
    end
    savefield = 'cohscal_ff_rho_rhoSE_base';
    
    % ###################### RUN
    ff = OUTSTRUCT.ffvals{i}(indsthis);
    cohscal = OUTSTRUCT.cohscal{i}(indsthis);
    
    Qfunc = @(x)corr(x(:,1), x(:,2), 'type', corrtype);
    X = [ff', cohscal'];
    rho_base = Qfunc(X);
    
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
    
    
    % ============== SAVE
    OUTSTRUCT.(savefield) = [OUTSTRUCT.(savefield); ...
        [rho_base rho_SE]];

    
    %% ============ WN corr
    if useonlybaseepoch==1
        indsthis = OUTSTRUCT.indsWN_epoch{i};
    else
        indsthis = OUTSTRUCT.indsWN{i};
    end
    savefield = 'cohscal_ff_rho_rhoSE_WN';
    
    % ###################### RUN
    ff = OUTSTRUCT.ffvals{i}(indsthis);
    cohscal = OUTSTRUCT.cohscal{i}(indsthis);
    
    Qfunc = @(x)corr(x(:,1), x(:,2), 'type', corrtype);
    X = [ff', cohscal'];
    rho_base = Qfunc(X);
    
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
    
    
    % ============== SAVE
    OUTSTRUCT.(savefield) = [OUTSTRUCT.(savefield); ...
        [rho_base rho_SE]];
    
    

    %% ============ All trials corr
    if useonlybaseepoch==1
        indsthis = [OUTSTRUCT.indsbase_epoch{i} OUTSTRUCT.indsWN_epoch{i}];
    else
        indsthis = [find(OUTSTRUCT.indsbase{i}) find(OUTSTRUCT.indsWN{i})];
    end
    savefield = 'cohscal_ff_rho_rhoSE_All';
    
    % ###################### RUN
    ff = OUTSTRUCT.ffvals{i}(indsthis);
    cohscal = OUTSTRUCT.cohscal{i}(indsthis);
    
    Qfunc = @(x)corr(x(:,1), x(:,2), 'type', corrtype);
    X = [ff', cohscal'];
    rho_base = Qfunc(X);
    
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
    
    
    % ============== SAVE
    OUTSTRUCT.(savefield) = [OUTSTRUCT.(savefield); ...
        [rho_base rho_SE]];

    %% ============ WN minus base, change in coh
    cohscal = OUTSTRUCT.cohscal{i};
    indsWN = OUTSTRUCT.indsWN_epoch{i};
    indsbase = OUTSTRUCT.indsbase_epoch{i};
    if cohdiff_usedprime==0
        cohscal_WnMinusBase = mean(cohscal(indsWN) - mean(cohscal(indsbase)));
    else
        asdfasdf;
    end
    
    % ============== SAVE
    OUTSTRUCT.cohscal_WN_minus_base = [OUTSTRUCT.cohscal_WN_minus_base; ...
        cohscal_WnMinusBase];
    
end