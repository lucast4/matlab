function [datbase, datWN, datWN_notminshuff, datbase_notminshuff, Xq] = ...
    lt_neural_POPLEARN_XCov_sub1(datmat_real, datmat_shuff, dosmooth, ...
    dosmooth_sigma, inds_base, inds_WN, xbins)

%% ==================== SMOOTH?
if dosmooth==1
    % SMOOTH PARAMS (gaussian)
    N = 2*dosmooth_sigma/0.001;
    win = gausswin(N);
    win = win./sum(win);
    % INTERP PARAMS.
    Xq = xbins(1):0.001:xbins(end);
    
    % 1) INTERPOLATE
    try
        tmp = interp1(xbins, datmat_real', Xq);
    catch err
        tmp = interp1(lags_sec, datmat_real', Xq);
    end
    
    datmat_real = tmp';
    % 2) SMOOTH
    for kkk=1:size(datmat_real,1)
        datmat_real(kkk,:) = conv(datmat_real(kkk, :), win, 'same');
    end
    
    % 1) INTERPOLATE
    tmp = interp1(xbins, datmat_shuff', Xq);
    datmat_shuff = tmp';
    % 2) SMOOTH
    for kkk=1:size(datmat_shuff,1)
        datmat_shuff(kkk,:) = conv(datmat_shuff(kkk, :), win, 'same');
    end
else
    Xq = [];
end

%% get xcov
datbase = mean(datmat_real(inds_base,:)) ...
    - mean(datmat_shuff(inds_base,:));

datWN = mean(datmat_real(inds_WN,:)) ...
    - mean(datmat_shuff(inds_WN,:));

datWN_notminshuff = mean(datmat_real(inds_WN,:));
datbase_notminshuff = mean(datmat_real(inds_base,:));

