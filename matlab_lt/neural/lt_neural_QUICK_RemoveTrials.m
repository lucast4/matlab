function badtrials = lt_neural_QUICK_RemoveTrials(bname, ename, sw, tvals, ...
    removebadtrialtype)
%% lt 2/5/19 - removes by hand bad trials, designated by time bnelow

% Enter array of trials (by datenums) and will tell you which ones to keep.
%^ NOTEL: can have multiple entries for a given bird expt pair. will
%incorporate all (i.e.give you the total badness).


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
    {{'pu69wh78', 'RALMANlearn1', 1, {'01Nov2017-2015', '01Nov2017-2120'}}}; ...
    {{'pu69wh78', 'RALMANlearn2', 1, {'05Nov2017-1130', '05Nov2017-1650'}}}; ...  05Nov2017-165006, ch21 goes bad after WN onset. % NOTE: looked thru represtentativbe of all song up to 1811 (i.e. switch to opposite direction). strange that there is island of songs that are OK, analyze that.
    {{'pu69wh78', 'RALMANlearn2', 1, {'05Nov2017-1733', '05Nov2017-2050'}}}; ... % ch21 goes bad after WN onset % 05Nov2017-173327
    {{'wh44wh39', 'RALMANlearn3', 1, {'21Mar2018-1759', '21Mar2018-2105'}}}; ... % 15 gets much worse, and others get worse too, but not as bad...
    {{'gr48bu5', 'RALMANLearn5', 1, {'23Jan2019-1724', '23Jan2019-2103'}}}; ... % then gets noisy at a(j). aj(j), and a(b), and j(b)...
    {{'gr48bu5', 'RALMANLearn6', 1, {'28Jan2019-2039', '28Jan2019-2100'}}}, ... % could keep, but gets a little bit noisier, and already showing good learning, so take earlier...
    ];
end

% ============== LFP NOISY (LARGE SHARED DEFLECTIONS)
if strcmp(removebadtrialtype, 'lfp') | strcmp(removebadtrialtype, 'all')
TrialsBad = [TrialsBad; ...
    {{'wh44wh39', 'RALMANlearn1', 1, {'19Feb2018-1755', '19Feb2018-2300'}}}; ...
    {{'wh72pk12', 'RALMANLearn3', 11, {'10Dec2018-2053', '10Dec2018-2310'}}}]; % actually ends at 1821, but just extend to end of day.
end

%% ========== ask whether the input matches any of the bad CHANS
badtrials = zeros(size(tvals)); % start with all good (all 0), then
for j=1:length(TrialsBad)
    tmp1 = strcmp(TrialsBad{j}{1}, bname);
    tmp2 = strcmp(TrialsBad{j}{2}, ename);
    tmp3 = TrialsBad{j}{3}==sw;
    
    if all([tmp1 tmp2 tmp3])
        disp(['BAD TRIALS: ' bname '-' ename '-sw' num2str(sw)]);
        
        % === figur eout which trials bad
        ton = datenum(TrialsBad{j}{4}{1}, 'ddmmmyyyy-HHMM');
        toff = datenum(TrialsBad{j}{4}{2}, 'ddmmmyyyy-HHMM');
        assert(toff>ton);
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

