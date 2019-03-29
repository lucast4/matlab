function chanbad = lt_neural_QUICK_RemoveBadChans(bname, ename, swnum, chan)
%% lt 3/21/19 - modified so doens't have to be just learning experiments.
% keep switch empty

%% lt 2/4/19 - tells you whether this is a bad channel.

% NOTE: chan can be array of channels. 
% NOTE: chan (within ChansBad) can also be array of channels)
%% NOTE DOWN BAD CHANS
% @# = NOIASE, base and WN
% @ = NOise, base
% # = NOISE, WN
% ? = potentially salvagable (i.e. only noisy epochs, could remove those)

%  THIS IS GOOD VERSION, use for general cases (e.g NGRAM)
ChansBad = {...,
    {'pu69wh78', 'RALMANlearn1', 1, 14}, ...  % # 
    {'pu69wh78', 'RALMANOvernightLearn1', [], 21}, ...  % 
    {'gr48bu5', 'RALMANLearn5', 1, 21}, ...  % # = 20 and 21 seem almost identical, so only keep 1. 20 seems better.
    {'wh44wh39', 'RALMANlearn2', [], 15}, ...  % # see lt_neural_QUICK_RemoveBadSyls
    {'wh44wh39', 'RALMANlearn4', [], 14}, ...  % # see lt_neural_QUICK_RemoveBadSyls
    {'wh72pk12', 'RALMANLearn3', [], 20}, ...  % # see lt_neural_QUICK_RemoveBadSyls - reason, not syl modulated (some song modualtion). interneuron?
    };



%% ========== ask whether the input matches any of the bad CHANS
chanbad = 0;
for j=1:length(ChansBad)
    tmp1 = strcmp(ChansBad{j}{1}, bname);
    tmp2 = strcmp(ChansBad{j}{2}, ename);
    tmp3 = ChansBad{j}{3}==swnum;
    tmp4 = any(ismember(ChansBad{j}{4}, chan));
    
    if all([tmp1 tmp2 tmp3 tmp4])
        chanbad=1;
        disp(['BAD CHAN: ' bname '-' ename '-sw' num2str(swnum) '-' num2str(chan)]);
        break
    end
end


