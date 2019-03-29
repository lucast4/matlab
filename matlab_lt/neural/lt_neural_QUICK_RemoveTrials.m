function badtrials = lt_neural_QUICK_RemoveTrials(bname, ename, sw, tvals, ...
    removebadtrialtype, chanthis)
%% lt 3/21/19 - modified so that can work even if not learnig experiment
% make the switch input empty.

%% lt 2/5/19 - removes by hand bad trials, designated by time bnelow

% Enter array of trials (by datenums) and will tell you which ones to keep.
%^ NOTEL: can have multiple entries for a given bird expt pair. will
%incorporate all (i.e.give you the total badness).

if ~exist('chanthis', 'var')
    chanthis = [];
end

if ~exist('removebadtrialtype', 'var')
    removebadtrialtype = 'spikes';
end

% if isempty(removebadtrialtype)
%     removebadtrialtype = 'spikes';
% end
%% enter bad trials here, separated by time.
% times are inclusive.

TrialsBad = {};

% ============== SPIKING BAND NOISY TRIALS
if strcmp(removebadtrialtype, 'spikes') | strcmp(removebadtrialtype, 'all')
    TrialsBad = [TrialsBad;
        {{'pu69wh78', 'RALMANOvernightLearn1', 1, {'10Nov2017-1156', '10Nov2017-2300'}, [9]}}; ... % RA channel is bad..
        {{'pu69wh78', 'RALMANlearn1', 1, {'01Nov2017-2015', '01Nov2017-2120'}, []}}; ...
        {{'pu69wh78', 'RALMANlearn2', 1, {'05Nov2017-1000', '05Nov2017-1002'}, []}}; ...  05Nov2017-165006, ch21 goes bad after WN onset. % NOTE: looked thru represtentativbe of all song up to 1811 (i.e. switch to opposite direction). strange that there is island of songs that are OK, analyze that.
        {{'pu69wh78', 'RALMANlearn2', 1, {'05Nov2017-1047', '05Nov2017-1047'}, []}}; ...  05Nov2017-165006, ch21 goes bad after WN onset. % NOTE: looked thru represtentativbe of all song up to 1811 (i.e. switch to opposite direction). strange that there is island of songs that are OK, analyze that.
        {{'pu69wh78', 'RALMANlearn2', 1, {'05Nov2017-1130', '05Nov2017-1650'}, []}}; ...  05Nov2017-165006, ch21 goes bad after WN onset. % NOTE: looked thru represtentativbe of all song up to 1811 (i.e. switch to opposite direction). strange that there is island of songs that are OK, analyze that.
        {{'pu69wh78', 'RALMANlearn2', 1, {'05Nov2017-1725', '05Nov2017-2050'}, []}}; ... % ch21 goes bad after WN onset % 05Nov2017-173327 [USED TO BE 1733]
        {{'wh44wh39', 'RALMANlearn1', 1, {'19Feb2018-1650', '19Feb2018-1704'}, []}}; ... % 15 gets much worse, and others get worse too, but not as bad...
        {{'wh44wh39', 'RALMANlearn3', 1, {'21Mar2018-1759', '21Mar2018-2105'}, [17 18 20 22]}}; ... % 15 gets much worse, and others get worse too, but not as bad...
        {{'gr48bu5', 'RALMANLearn5', 1, {'23Jan2019-1724', '23Jan2019-2103'}, []}}; ... % then gets noisy at a(j). aj(j), and a(b), and j(b)...
        {{'gr48bu5', 'RALMANLearn6', 1, {'28Jan2019-2039', '28Jan2019-2100'}, []}}; ... % could keep, but gets a little bit noisier, and already showing good learning, so take earlier...
        {{'wh72pk12', 'RALMANLearn7', 1, {'11Jan2019-2040', '11Jan2019-2300'}, []}}; ... % gradual fall off, see noise becuase spikes wwaker. this is rleatively conservative I think. could prb go to aroudn2140;
        {{'br92br54', 'LMANlearn6', [], {'16Jun2017-1146', '16Jun2017-1220'}, []}}; ... % see lt_neural_QUICK_RemoveBadSyls
        {{'or74bk35', 'LMANneural2', [], {'12May2017-1444', '12May2017-2133'}, []}}; ... % see lt_neural_QUICK_RemoveBadSyls
        ];
end

% TO ADD:
% {{'gr48bu5', 'RALMANLearn4', 1, {'20Jan2019-2000', '20Jan2019-2300'}, []}}; ... % ADDED NEW, gradual dropoff, but seems more noisy
% {{'gr48bu5', 'RALMANLearn5', 1, {'23Jan2019-1214', '23Jan2019-2300'}, []}}; ... % ADDED NEW - remove this experiment since (1) cannot tell if really spikes in RA, and (2) similar across RA sites. Baseline seems ok. Keep, since seems like is real spikes and is OK

% last ind (e.g. [17 18 20 22]) are channels to ignore for this - i.e. if
% 17, then wil not remove those trials.


% ============== LFP NOISY (LARGE SHARED DEFLECTIONS)
if strcmp(removebadtrialtype, 'lfp') | strcmp(removebadtrialtype, 'all')
    TrialsBad = [TrialsBad; ...
        {{'wh44wh39', 'RALMANlearn1', 1, {'19Feb2018-1755', '19Feb2018-2300'}, []}}; ...
        {{'wh72pk12', 'RALMANLearn3', 11, {'10Dec2018-2053', '10Dec2018-2310'}, []}}]; % actually ends at 1821, but just extend to end of day.
end

%% ========== ask whether the input matches any of the bad CHANS
badtrials = zeros(size(tvals)); % start with all good (all 0), then
for j=1:length(TrialsBad)
    tmp1 = strcmp(TrialsBad{j}{1}, bname);
    tmp2 = strcmp(TrialsBad{j}{2}, ename);
    tmp3 = TrialsBad{j}{3}==sw;
    
    if ~isempty(chanthis)
        if any(ismember(TrialsBad{j}{5}, chanthis))
            % then ignore - this channel is exception, so not bad trials.
            continue
        end
    end
    
    if all([tmp1 tmp2 tmp3])
        disp(['BAD TRIALS: ' bname '-' ename '-sw' num2str(sw)]);
        
        % === figur eout which trials bad
        ton = datenum(TrialsBad{j}{4}{1}, 'ddmmmyyyy-HHMM');
        toff = datenum(TrialsBad{j}{4}{2}, 'ddmmmyyyy-HHMM');
        assert(toff>=ton);
        % --- which trials are within these bounds?
        tmp = tvals>=ton & tvals<=toff;
        badtrials = badtrials | tmp; % keep a running trakc of all bad triasl...
    end
end

badtrials = logical(badtrials);

% inds1 = cellfun(@(x)strcmp(x{1}, bname), TrialsBad);
% inds2 = cellfun(@(x)strcmp(x{2}, ename), TrialsBad);
% indthis = find(inds1 & inds2);
%
% if ~any(indthis)
%     badtrials = zeros(size(tvals));
%     return
% elseif sum(indthis)>1
%     disp('WHY multile? [check entries not duplicated...');
%     pause
% else
%    % ===== figure out which trials to keep
%    TrialsBad
% end
%

