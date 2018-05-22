function lt_neural_v2_CTXT_BRANCH_DatVsShuff(analyfname, Niter, TimeWindows, ...
    dotransform)
%% 5/11/18 - modified to allow for both alignment by onset and offset, detects automatically from params

% if align by onset: then will take from sylonset+TimeWindows(1) to
% syloffset (i.e. 25th percentile of syl dur) + TimeWindows(2)

% if align by offset: then will take from syloffset - syldur(25th
% percentile) + TimeWindows(1) to offset + TimeWindows(2)
% so if TimeWindows = [-0.01 -0.01] then will get from 10ms before syl
% onset (infered as syl offset minus 25th percentile of syl dur) to 10ms
% before syl offets.

%% lt 11/28/17 - output saved in a directory (as individual .mat files)
% PREVIOUSLY saved as a field in the CLASSES structure, but it is too
% large, so do this instead.

%% lt 10/26/17 - for premotor window, calculates decoding accuracy.
% TO implement: sliding time window, assuming a certain delay between LMAN
% and RA


%% load branch
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

load([savedir '/CLASSESv2_' analyfname '.mat']);
load([savedir '/SUMMARYv2_' analyfname '.mat']);
try
    load([savedir '/PARAMSv2_' analyfname '.mat']);
catch err
    
end

% assert(length(ALLBRANCH.alignpos)==1, 'note: using old version of branch, since has >1 alignpos in one struct ..., code wil only do alignpos(1)');

%% params

% TimeWindows = [-0.05 -0.05]; % [-0.05 -0.05] means window from 50ms pre onset to 50ms pre offset (each row is separate analysis)
% TimeWindows = [-0.05 -0.05;-0.1 -0.1]; % [-0.05 -0.05] means window from 50ms pre onset to 50ms pre offset (each row is separate analysis)
Binsize = [0.005]; % 5 ms.
numtimebins = size(TimeWindows,1);

%
prms.ClassSlide.frbinsize =Binsize;
prms.ClassSlide.Nmin = 7;

%%- decoding
CVkfoldnum = min([prms.ClassSlide.Nmin, 5]); % seems comparable to LOO at 8fold.
rebalance =1;
imbalance_thr = 0.7;
beta = 0.9;

% Niter = 100; % num shuffles of neg control
decodestat = 'F1';


%%

currdir = pwd;
try
    cd([savedir '/' analyfname '/SHUFFDECODE/']);
    cd(currdir)
catch err
    mkdir([savedir '/' analyfname '/SHUFFDECODE/']);
end

%% RUN

numbirds = length(CLASSES.birds);

for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    
    for ii=1:numneurons
        
        numbranch = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        for iii=1:numbranch
            
            disp(['brd' num2str(i) '-n' num2str(ii) '-br' num2str(iii)]);
            
            datstruct = CLASSES.birds(i).neurons(ii).branchnum(iii);
            
            
            % ======= to get offset of window, take deviation
            % from 25th percentile (short side) offset of all contexts in this branch
            numclasses = length(datstruct.SEGEXTRACT.classnum);
            alldurs = [];
            for cc = 1:numclasses
                if isempty(datstruct.SEGEXTRACT.classnum(cc).SegmentsExtract)
                    continue
                end
                tmp = [datstruct.SEGEXTRACT.classnum(cc).SegmentsExtract.Dur_syl];
                alldurs = [alldurs tmp];
            end
            
            % ============= what to use as syl dur depends on whether
            % aligning to onset or offest
            syldur = prctile(alldurs, 25);
            
            TimeWindows_relonset = TimeWindows;
            if prms.alignOnset==1
                % then straightforward - get time relative to onset
            TimeWindows_relonset(:,2) = TimeWindows_relonset(:,2)+syldur; % converted... can easily modify this iof wanted
            elseif prms.alignOnset==0
                % get times relative to syllable offset ---
                TimeWindows_relonset(1) = TimeWindows_relonset(1) - syldur;
            end
            
            
            % ################ go thru all time bins
            for tt = 1:numtimebins
                
                % ================ check if already done - if so, skip
                savefname = [savedir '/' analyfname '/SHUFFDECODE/bird' num2str(i) '_neur' num2str(ii) ...
                    '_branch' num2str(iii) '_tbin' num2str(tt) '.mat'];
                if exist(savefname, 'file')
                    disp(['SKIP - already exist: ' savefname]);
                    continue
                end
                
                % ------ SOME PARAMS
                prms.classtmp.frtimewindow = TimeWindows_relonset(tt,:); % on and off, relative to syl onset
                
                % ================= GET DATA IN FORMAT FOR CLASSIFICATION
                SEGEXTRACT = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT;
                if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                    clustnum = [];
                else
                    clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                end
                
                [Xall, ~, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum, ...
                    dotransform);
                
                if length(CtxtClasses)<2
                    continue
                end
                
                % ============================ get decode of actual data
                NNN = 3;
                ConfMatAll = cell(1,NNN);
                for nn = 1:NNN
                    % do 3 times and take mean
                    
                    [Ypredicted, ConfMat] ...
                        = lt_neural_v2_QUICK_classify(Xall, Y, 'glmnet', ...
                        rebalance, imbalance_thr, beta, CVkfoldnum);
                    
                    ConfMatAll{nn} = int8(ConfMat);
                end
                
                
                
                % =========================== get shuffled decodes for this
                % time bin
                ConfMatAll_NEG = cell(1,Niter);
                for nn=1:Niter
                    indtmp = randperm(size(Y,1));
                    Yperm = Y(indtmp);
                    
                    % --- classify
                    [Ypredicted, ConfMat] ...
                        = lt_neural_v2_QUICK_classify(Xall, Yperm, 'glmnet', ...
                        rebalance, imbalance_thr, beta, CVkfoldnum);
                    
                    assert(~isempty(Ypredicted), 'whye mpty?');
                    
                    ConfMatAll_NEG{nn} = int8(ConfMat);
                end
                
                % ================= OUTPUT
                if (0)
                    % old version
                    CLASSES.birds(i).neurons(ii).branchnum(iii).SHUFFDECODE.timebin(tt).window_relonset = prms.classtmp.frtimewindow;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).SHUFFDECODE.timebin(tt).ConfMatAll_DAT = ConfMatAll;
                    CLASSES.birds(i).neurons(ii).branchnum(iii).SHUFFDECODE.timebin(tt).ConfMatAll_NEG = ConfMatAll_NEG;
                    
                    
                end
                
                %% =================== compare data to distribution
                % -- get mean F1 for dat
                tmp = [];
                for j=1:length(ConfMatAll)
                    cmat = ConfMatAll{j};
                    sts = lt_neural_ConfMatStats(cmat);
                    tmp = [tmp sts.(decodestat)];
                end
                decode_dat = mean(tmp);
                
                % -- compare to negative control distr
                decode_neg = nan(1,length(ConfMatAll_NEG));
                for j=1:length(ConfMatAll_NEG)
                    cmat = ConfMatAll_NEG{j};
                    
                    sts = lt_neural_ConfMatStats(cmat);
                    decode_neg(j) = sts.(decodestat);
                end
                
                if (0) % visualize
                    lt_figure; hold on;
                    lt_plot_histogram(decode_neg);
                    line([decode_dat decode_dat], ylim, 'Color', 'r');
                end
                
                % --- prob of dat in null
                Pdat = sum(decode_neg>=decode_dat)./length(decode_neg);
                
                
                % ==================== OUTPUT
                if (0)
                    % -- old version
                    CLASSES.birds(i).neurons(ii).branchnum(iii).SHUFFDECODE.timebin(tt).Pdat = Pdat;
                end
                
                % ================= NEW VERSION - SAVES EACH OUTPUT
                decodestruct = struct;
                decodestruct.window_relonset = prms.classtmp.frtimewindow;
                decodestruct.ConfMatAll_DAT = ConfMatAll;
                decodestruct.ConfMatAll_NEG = ConfMatAll_NEG;
                decodestruct.Pdat = Pdat;
                
                save(savefname, 'decodestruct');
                
                
            end
        end
    end
end

%% ========== save params

savename_par = [savedir '/' analyfname '/SHUFFDECODE/Params.mat'];
save(savename_par, 'TimeWindows');

%% ======= save classes (overwrite old struct)
if (0)
    CLASSES.SHUFFDECODEpar.TimeWindows_relOnsetOffset =TimeWindows;
    save([savedir '/CLASSESv2_' analyfname '.mat'], 'CLASSES');
end

%% ================== debug, to convert from older bersion (saving in struct) to new version (saving .mat)
if (0)
    lt_neural_trash;
end
end


function [Xall, xtimesall, Y, CtxtClasses] = fn_extractClassDat(SEGEXTRACT, prms, clustnum, ...
    dotransform)

frtimewindow = prms.classtmp.frtimewindow; % on and off, relative to syl onset
frbinsize = prms.ClassSlide.frbinsize;
Nmin = prms.ClassSlide.Nmin;

numclasses = length(SEGEXTRACT.classnum);

% =================== COLLECT DATA TO CLASSIFY
Xall = []; % trials x bins (FR vectors)
xtimesall = []; % 1 x bins (time stamps)
Y = []; % context indicator
CtxtClasses = {};
contextcounter = 0;


% ------ methiod2
for j=1:numclasses
    
    sylname = SEGEXTRACT.classnum(j).regexpstr;
    
    % --- EXTRACT DATA
    segextract = SEGEXTRACT.classnum(j).SegmentsExtract;
    
    % -- extract FR
    segextract = lt_neural_SmoothFR(segextract, clustnum);
    
    
    if ~isfield(segextract, 'spk_Times')
        % then no data
        continue
    end
    
    if length(segextract) < Nmin
        % not neough data
        continue
    end
    
    
    % ---- EXTRAC FR VECTOR (within desired time window)
    xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
    
    % --------- min val depends on whether aligned to onset or offset
    xminval = prms.motifpredur+frtimewindow(1)+0.0001;
    xmaxval = prms.motifpredur+frtimewindow(2)+0.0001;
    indsFR = xbin>(xminval) & xbin<=(xmaxval); % add/minus at ends because somtimes
    
    X = [segextract.FRsmooth_rate_CommonTrialDur];
    X = X(indsFR, :);
    X = X';
    
    xtimes = xbin(indsFR);
    xtimes = xtimes';
    
    
    % ------------ reduce Dim of FR vectors (by binning in
    % time)
    TrimDown = 1;
    [X, xtimes] = lt_neural_v2_QUICK_binFR(X, xtimes, frbinsize, TrimDown);
    
    % ==================== do square root transform?
    if dotransform==1
        X = sqrt(X);
    end
    
    
    % ======================== COLLECT ACROSS ALL CLASSES
    contextcounter = contextcounter+1;
    
    CtxtClasses = [CtxtClasses sylname];
    
    Xall = [Xall; X];
    xtimesall = [xtimesall; xtimes];
    Y = [Y; contextcounter*ones(size(X,1),1)];
    
end

% --- convert Y to categorical array
if version('-release')=='2013a'
    Y = nominal(Y);
else
    Y = categorical(Y);
end
end