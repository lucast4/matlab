function [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds, ...
    DecodeDAT, DecodeShuff, Decode_Z] = lt_neural_NGRAMS_QUICK_FRdiff(frmat1, frmat2, downfactor, windowx, ...
    Nmin, nshufftmp, Params, DoDecode)
%% added decode 
% DoDecode = 1; % to get decoder performances; IMPORTNAT: IF 1, then overwrites 
% FRdiff values.

%% lt 5/4/18 - calculate diff in FR, given frmats

% e.g. fr mats are
% frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat; % timebin x trial
% [RAW]
% frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
% windowx = motifpredur+Params.window_prem; % i.e all datapoints >=0. windwo to take FR
% Nmin = Params.Nmin; i.e. downsample will make sure to not go below this.
% i.e. downsample is done on trials in excess of this.

% ---------------- 2) DOWNSAMPLE
% downfactor = 0.5; % to downsample, while still at least equal to Nmin
% leave empty to not downsample

% ---------------
% DoDecode = 1, then does logistic regression
%% ============
% since will get zscore, shuff should be at least 2
nshufftmp = max([2 nshufftmp]);


%% ========== prepare fr mat

% ---------------- original sample sizes;
N1 = size(frmat1, 2);
N2 = size(frmat2, 2);

% ---------------- 2) DOWNSAMPLE
if downfactor <1
    N1_new = round(downfactor*(N1 - Nmin)) + Nmin; % downfactor*[excess from Nmin] + [Nmin]
    N2_new = round(downfactor*(N2 - Nmin)) + Nmin; % downfactor*[excess from Nmin] + [Nmin]
    
    % ================= replace things
    frmat1 = frmat1(:, randperm(N1, N1_new));
    frmat2 = frmat2(:, randperm(N2, N2_new));
    N1 = N1_new;
    N2 = N2_new;
end

% =========== Concatenate into FRmat
FRmat = [frmat1 frmat2];
x = lt_neural_QUICK_XfromFRmat(frmat1);
TrialInds = {1:N1, N1+1:N1+N2};

% -------------- pare down to just premotor windwo
FRmat = FRmat(x>=windowx(1) & x<=windowx(2),:);

% ------------------- TAKE SQUARE ROOT TRANSFORM
if Params.doSqrtTransform==1
    FRmat = sqrt(FRmat);
end

% ============= OUTPUT
Nboth = [N1; N2];

%% ############################## FR DIFF (DAT)
fr1 = mean(FRmat(:, TrialInds{1}),2);
fr2 = mean(FRmat(:, TrialInds{2}),2);
var1 = var(FRmat(:, TrialInds{1}), 0, 2);
var2 = var(FRmat(:, TrialInds{2}), 0, 2);

if Params.use_dPrime==1
    dpr = (fr1-fr2)./sqrt(0.5*(var1 + var2));
    FRdiffDAT = mean(abs(dpr));
else
    % then do regular FR diff
    FRdiffDAT = mean(abs(fr1-fr2));
end

%% ############################ CONTROLS
FRdiffShuff = [];

nsamps = size(FRmat,2);
for ss = 1:nshufftmp
    
    indshuff = randperm(nsamps);
    FRmatSHUFF = FRmat(:,indshuff);
    
    % ================= calculate
    fr1 = mean(FRmatSHUFF(:, TrialInds{1}),2);
    fr2 = mean(FRmatSHUFF(:, TrialInds{2}),2);
    var1 = var(FRmatSHUFF(:, TrialInds{1}), 0, 2);
    var2 = var(FRmatSHUFF(:, TrialInds{2}), 0, 2);
    
    if Params.use_dPrime==1
        dpr = (fr1-fr2)./sqrt(0.5*(var1 + var2));
        frdiff = mean(abs(dpr));
    else
        % then do regular FR diff
        frdiff = mean(abs(fr1-fr2));
    end
    
    FRdiffShuff = [FRdiffShuff; frdiff];
end

%% ########################## zscore data relative to shuff
FRdiff_Z = (FRdiffDAT - mean(FRdiffShuff))./std(FRdiffShuff);

%% ######################## DECODE
if DoDecode==1
    
    % ----- 1) BIN FR MATRIX
    [Xall, xtimesall, Y] = ...
        lt_neural_NGRAMS_DecodePrep(frmat1, frmat2, windowx, ...
        Params);
    
    % ------ 2) DECODER
    rebalance =1;
    imbalance_thr = 0.7;
    beta = 0.9;
    CVkfoldnum = 4;
    [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
        = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
        rebalance, imbalance_thr, beta, CVkfoldnum);
    
    % ---------- OUTPUT
    sts = lt_neural_ConfMatStats(ConfMat);
    DecodeDAT = sts.F1;
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHUFFLE
    FRmatRaw = [frmat1 frmat2];
    nsamps = size(FRmatRaw,2);
    DecodeShuff = [];
    for ss = 1:nshufftmp
        
        indshuff = randperm(nsamps);
        FRmatSHUFF = FRmatRaw(:, indshuff);
        
        % ================= calculate
        frmat1SHUFF = FRmatSHUFF(:, TrialInds{1});
        frmat2SHUFF = FRmatSHUFF(:, TrialInds{2});
        
        
        % ################################## DECODE
        [Xall, xtimesall, Y] = ...
            lt_neural_NGRAMS_DecodePrep(frmat1SHUFF, frmat2SHUFF, windowx, ...
            Params);
        
        % ------ 2) DECODER
        rebalance =1;
        imbalance_thr = 0.7;
        beta = 0.9;
        CVkfoldnum = 4;
        [Ypredicted, ConfMat, accuracy, sensitivity_mean, PiYActual] ...
            = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
            rebalance, imbalance_thr, beta, CVkfoldnum);
        
        % #############################################
        sts = lt_neural_ConfMatStats(ConfMat);
        DecodeShuff = [DecodeShuff; sts.F1];
    end
    
    Decode_Z = (DecodeDAT - mean(DecodeShuff))./std(DecodeShuff);
    
    
    
    % ====================== replace output with decode
    FRdiffDAT = DecodeDAT;
    FRdiffShuff = DecodeShuff;
    FRdiff_Z = Decode_Z;
    
end




