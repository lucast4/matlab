function sylbad = lt_neural_QUICK_RemoveBadSyl(bname, ename, syltoken, ...
    typesToRemove, chanthis)
%% lt 3/21/19 - modified lt_neural_QUICK_LearnRemoveBadSyl to be general across experiments
% 3/21/19 - UP TO DATE - all non-learning data (i.e. if learning, then is
% preceding first switch).

% METHOD: by eye looked thru bunch of songs and if see consistent noise
% then flag for removal.
% NOTE: if a channel is good, then give it an exception.

% NOTE: not based at all on lt_neural_QUICK_LearnRemoveBadSyl, since the latter looks thru both WN and baseline.

%% NOTES (3/21/19) - these are notes taken as I was lookgin thruogh. 
% use the following to make this list - i.e. for each bird plot each unqiue
% directory + channels that are in that directory:

% longversion = 0; % then shows each extracted neuron
% lt_neural_SummaryStruct_ShowDirChans(SummaryStruct, longversion);


% /bluejay4/lucas/birds/bk7/NEURAL/080516_Neural_Audio  ---  14  19  23
% /bluejay4/lucas/birds/bk7/NEURAL/080616_Neural_Audio  ---  9  14
% /bluejay4/lucas/birds/bk7/NEURAL/080916_NeuralAudio_  ---  9  14

% /bluejay4/lucas/birds/bk7/NEURAL/081616_neuralAudio  ---  18 - slight
% signal loss

% /bluejay4/lucas/birds/bk7/NEURAL/081916_neuralaudio  ---  14  18 slight signal loss
% /bluejay4/lucas/birds/bk7/NEURAL/082216_neuralaudio  ---  23 slight signal loss
% /bluejay4/lucas/birds/bk7/NEURAL/100516_LearnLMAN1  ---  10  18  23
% /bluejay4/lucas/birds/bk7/NEURAL/100616_LearnLMAN1  ---  10 n1, r, remove? sometimes noise
% /bluejay4/lucas/birds/bk7/NEURAL/101716_LearnLMAN2  ---  10  14, both noisy. n1, r, (j)k, e, u, remove

% /bluejay5/lucas/birds/bu77wh13/NEURAL/020317_LMANneural1  ---  9 GOOD
% /bluejay5/lucas/birds/bu77wh13/NEURAL/020417_LMANneural1  ---  9  18
% /bluejay5/lucas/birds/bu77wh13/NEURAL/020517_LMANlearn1  ---  9
% /bluejay5/lucas/birds/bu77wh13/NEURAL/020517_LMANneural1  ---  9 g(k) bad
% /bluejay5/lucas/birds/br92br54/NEURAL/042017_LMANlearn2  ---  14  18  23
% /bluejay5/lucas/birds/br92br54/NEURAL/042317_LMANlearn3  ---  14  18  23
% /bluejay5/lucas/birds/br92br54/NEURAL/042717_LMANlearn4  ---  14  18  23 
% /bluejay5/lucas/birds/br92br54/NEURAL/042917_LMANlearn4  ---  14  18  23
% /bluejay5/lucas/birds/br92br54/NEURAL/052417_LMANlearn5  ---  14 s, d, y, m, all sometimes with hi hz noise, but not always, so seems fine.
% /bluejay5/lucas/birds/br92br54/NEURAL/061617_LMANlearn6  ---  14  23 remove after 1146. remove s(d), c(s), s(m), s(y)
% /bluejay5/lucas/birds/or74bk35/NEURAL/051217_LMANneural2  ---  8   9  17  20 REMOVE : too much shared noise

% /bluejay5/lucas/birds/pu69wh78/NEURAL/103017_RAlearn1  ---  17  18  21 OK
% /bluejay5/lucas/birds/pu69wh78/NEURAL/103117_LMANsearch  ---  9  14 (a1 remove some)
% /bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1  ---  9  14  17  18  21 OK
% /bluejay5/lucas/birds/pu69wh78/NEURAL/110517_RALMANlearn2  ---  9  21 OK
% /bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1  --- 9  14  21 K


% /bluejay5/lucas/birds/wh44wh39/NEURAL/021518_Rec3  ---  15  20  21  22 OK, remove (c1) for 1442 and after
% /bluejay5/lucas/birds/wh44wh39/NEURAL/021918_RALMANlearn1  ---  15  17  18  20  21  22 OK

% /bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2  ---  14  15  17  18  20  21 remove chan 15
% /bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3  ---  14  15  17  18  20  22 GOOD

% /bluejay5/lucas/birds/wh44wh39/NEURAL/041718_RALMANlearn4  ---  14  15  17  21  22 remove 14. 15: remove c1, (m)a, m(a), (j)j
% /bluejay0/bluejay2/lucas/birds/wh72pk12/NEURAL/010419_RALMANLearn6  ---  21
% /bluejay0/bluejay2/lucas/birds/wh72pk12/NEURAL/011119_RALMANLearn7  ---  21
% /bluejay0/bluejay2/lucas/birds/wh72pk12/NEURAL/011219_RALMANLearn7  ---  21
% /bluejay0/bluejay2/lucas/birds/wh72pk12/NEURAL/011619_RALMANLearn8  ---  21
% /bluejay0/bluejay2/lucas/birds/wh72pk12/NEURAL/113018_RALMANLearn3  ---  8  14  15  20  21 remove (j)kl. maybe j(r)b


% /bluejay0/bluejay2/lucas/birds/gr48bu5/NEURAL/011519_RALMANLearn2  ---  8  15  18  20  21  23 remove d(a) and h(a)
% /bluejay0/bluejay2/lucas/birds/gr48bu5/NEURAL/011719_RALMANLearn3  ---  15  20  21 remove d(a) and h(a)
% /bluejay0/bluejay2/lucas/birds/gr48bu5/NEURAL/012019_RALMANLearn4  ---  15  21 remove d(a) and h(a) but could keep, not too bad
% /bluejay0/bluejay2/lucas/birds/gr48bu5/NEURAL/012319_RALMANLearn5  ---  15  20  21 remove all a
% /bluejay0/bluejay2/lucas/birds/gr48bu5/NEURAL/012819_RALMANLearn6  ---  8  15  20  23 remove all a

%% 
% chanthis is current channel. can exclude certain channels below.
if ~exist('chanthis', 'var')
    chanthis = [];
end

if ~exist('typesToRemove', 'var')
   typesToRemove = {'wn', 'noise'}; % by default remove all 
end

%% syllabels that follow the target (on same motif up to 2+ syls) should be excluded

% @# = NOIASE, base and WN
% @ = NOise, base
% # = NOISE, WN
% ? = potentially salvagable (i.e. only noisy epochs, could remove those)

% e.g. {'pu69wh78', 'RALMANOvernightLearn1', 1, 'aabh(h)', 'WN', [9 21]}, where
% "wn" or "noise" at end means that it is removed becuase it follows WN
% dring elarning or becuase of artifact (noise). if empty, then assues it
% is because follows "wn" [this is actually the case, since I did nto start
% including artifact cases until later]
% NOTE: "wn" could mean that is also has noise ...
% [9 21] refers to channels to make exception for. don't put anything to
% apply to all channels.
% entry like this ('0900-1441') means is good times (i.e. do not apply)

SylsBad = {...,
    {'bk7', 'LearnLMAN1', 'r', '', []}, ...
    {'bk7', 'LearnLMAN2', 'r', '', []}, ...
    {'bk7', 'LearnLMAN2', '(j)k', '', []}, ...
    {'bk7', 'LearnLMAN2', 'e', '', []}, ...
    {'bk7', 'LearnLMAN2', 'u', '', []}, ...
    {'bu77wh13', 'LMANneural1', 'g(k)', '', []}, ...
    {'br92br54', 'LMANlearn5', 's', '', []}, ...
    {'br92br54', 'LMANlearn5', 'd', '', []}, ...
    {'br92br54', 'LMANlearn5', 'y', '', []}, ...
    {'br92br54', 'LMANlearn5', 'm', '', []}, ...
    {'br92br54', 'LMANlearn6', 's(d)', '', []}, ...
    {'br92br54', 'LMANlearn6', 'c(s)', '', []}, ...
    {'br92br54', 'LMANlearn6', 's(m)', '', []}, ...
    {'br92br54', 'LMANlearn6', 's(y)', '', []}, ...
    {'pu69wh78', 'LMANsearch', '(a)a', '', []}, ...
    {'pu69wh78', 'RALMANlearn1', '(a)a', '', []}, ...
    {'pu69wh78', 'RALMANlearn2', '(a)a', '', []}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', '(a)a', '', []}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 'a(b)', '', []}, ...
    {'wh44wh39', 'Rec3', '(c)c', '', [], '0900-1441'}, ...
    {'wh44wh39', 'RALMANlearn3', 'm(j)', '', [14 17 18 20 22]}, ...
    {'wh44wh39', 'RALMANlearn4', '(c)c', '', [14 17 21 22]}, ...
    {'wh44wh39', 'RALMANlearn4', '(m)a', '', [14 17 21 22]}, ...
    {'wh44wh39', 'RALMANlearn4', 'm(a)', '', [14 17 21 22]}, ...
    {'wh44wh39', 'RALMANlearn4', '(j)j', '', [14 17 21 22]}, ...
    {'wh72pk12', 'RALMANLearn3', '(j)kl', ''}, ...
    {'gr48bu5', 'RALMANLearn2', 'd(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn2', 'h(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn3', 'd(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn3', 'h(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn4', 'd(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn4', 'h(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn5', 'd(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn5', 'h(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn6', 'd(a)', ''}, ...
    {'gr48bu5', 'RALMANLearn6', 'h(a)', ''}, ...
    };


%     {'wh72pk12', 'RALMANLearn6', 1, '(j)io', 'wn'}, ... % ************************************

%% ========== fill in whether is "noise" or "wn"
% === if empty, then is "wn"
for i=1:length(SylsBad)
   if length(SylsBad{i})==3
       % then add on the end
      SylsBad{i} = [SylsBad{i} 'noise'];
   elseif isempty(SylsBad{i}{4})
       SylsBad{i}{4} = 'noise';
   end
end


%% ========== ask whether the input matches any of the bad syls
sylbad = 0;
for j=1:length(SylsBad)
    tmp1 = strcmp(SylsBad{j}{1}, bname);
    tmp2 = strcmp(SylsBad{j}{2}, ename);
    tmp3 = strcmp(SylsBad{j}{3}, syltoken);
    
    % -- if this is chan to exclude, then abort
    if ~isempty(chanthis) & length(SylsBad{j})>4
        
        if any(ismember(SylsBad{j}{5}, chanthis))
        % then this channel should be excludede, i.e. this syl cannot be
        % bad for this chan
            continue
        end
    end
    
    if all([tmp1 tmp2 tmp3]) & any(ismember(typesToRemove, SylsBad{j}{4}))
        sylbad=1;
        disp(['BAD SYL: ' bname '-' ename '-' syltoken]);
        break
    end
end

