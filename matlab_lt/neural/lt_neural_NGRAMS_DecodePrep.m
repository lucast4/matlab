function [Xall, xtimesall, Y] = ...
                    lt_neural_NGRAMS_DecodePrep(frmat1, frmat2, windowx, ...
                    Params)
% 
% frmat1, frmat2: bin(0.001s) x trial

%% =================== COLLECT DATA TO CLASSIFY

Xall = []; % trials x bins (FR vectors)
xtimesall = []; % 1 x bins (time stamps)
Y = []; % context indicator
contextcounter = 0;

frbinsize = 0.005;


FRcell = {};
FRcell{1} = frmat1;
FRcell{2} = frmat2;

% ------ methiod2
for k=1:2
    
    xbin = lt_neural_QUICK_XfromFRmat(FRcell{k});
    
    % ---- EXTRAC FR VECTOR (within desired time window)
    indsFR = xbin>(windowx(1)+0.0001) ...
        & xbin<=(windowx(2)+0.0001); % add/minus at ends because somtimes
    
    X = FRcell{k};
    X = X(indsFR, :);
    X = X';
    
    xtimes = xbin(indsFR);
    
    
    % ------------ reduce Dim of FR vectors (by binning in
    % time)
    TrimDown = 1;
    [X, xtimes] = lt_neural_v2_QUICK_binFR(X, xtimes, frbinsize, TrimDown);
    
    
    % ==================== do square root transform?
    if Params.doSqrtTransform==1
        X = sqrt(X);
    end
    
    
    % ======================== COLLECT ACROSS ALL CLASSES
    contextcounter = contextcounter+1;
    
    Xall = [Xall; X];
    xtimesall = [xtimesall; xtimes];
    Y = [Y; contextcounter*ones(size(X,1),1)];
    
end

%% --- convert Y to categorical array
if version('-release')=='2013a'
    Y = nominal(Y);
else
    Y = categorical(Y);
end