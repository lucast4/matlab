function sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, swnum, ...
    syltoken, typesToRemove, chanthis)
%% 
% chanthis is current channel. can exclude certain channels below.
if ~exist('chanthis', 'var')
    chanthis = [];
end

if ~exist('typesToRemove', 'var')
   typesToRemove = {'wn', 'noise'}; % by default remove all 
end

%% tells you whether you should keep a given syl for a given switch, learing
% NOTE: have looked closely at all expts so far for wh44 and pu69


%%  TOOL TO LIST ALL SYLS THAT EXIST FOR A GIVEN SWITCH
if (0)
    for i=1:length(SwitchCohStruct.bird)
        bname = SwitchStruct.bird(i).birdname;
        for ii=1:length(SwitchCohStruct.bird(i).exptnum)
            ename = SwitchStruct.bird(i).exptnum(ii).exptname;
            
            for iii=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
                if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum)
                    continue
                end
                % ============================== get list of motifs.
                disp(' ======================================================= ');
                disp([bname ' -- ' ename  ' -- sw' num2str(iii)]);
                disp({SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum.motifname});
                
            end
        end
    end
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


SylsBad = {...,
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'aabh(h)', 'wn'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 1, '(a)ab', 'noise', [9]}, ... % @#
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'a(a)b', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'aa(b)', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'j(j)b', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'jj(b)', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'h(g)', 'wn', [9 14]}, ... % @#
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'jjbhh(g)', 'noise'}, ... % #
    {'pu69wh78', 'RALMANOvernightLearn1', 6, 'aab(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 6, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 6, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 6, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 7, 'aab(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 7, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 7, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 7, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 8, 'aab(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 8, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 8, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 8, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 9, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 9, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 10, 'aab(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 10, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 10, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 10, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 11, 'aab(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 11, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 11, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 11, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 12, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 12, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANlearn1', 1, '(a)ab', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANlearn2', 1, '(a)ab', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANlearn2', 1, 'a(a)b', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANlearn2', 1, 'j(j)bhh', 'noise'}, ... % @#
    {'pu69wh78', 'RALMANlearn2', 1, '(j)jbhh', 'noise'}, ... % #
    {'pu69wh78', 'RALMANlearn2', 1, 'jj(b)hh', 'noise'}, ... % #
    {'wh44wh39', 'RALMANlearn1', 1, 'dk(c)c', 'noise'}, ... % # [note: is actually not bad, except for ch17. others are sort of bad at around 1704+ (sometimes clack). Easiest is just to remove all c]
    {'wh44wh39', 'RALMANlearn1', 2, 'cb(b)'}, ...
    {'wh44wh39', 'RALMANlearn1', 7, 'cb(b)'}, ...
    {'wh44wh39', 'RALMANlearn1', 8, 'cb(b)'}, ...
    {'wh44wh39', 'RALMANlearn1', 9, 'cb(b)'}, ...
    {'wh44wh39', 'RALMANlearn2', 1, 'dk(c)c', 'noise'}, ... % @# all
    {'wh44wh39', 'RALMANlearn2', 1, '(m)d', 'noise', [18 20 21 22]}, ... % # 
    {'wh44wh39', 'RALMANlearn2', 1, '(n)hh', 'wn'}, ... % @# 
    {'wh44wh39', 'RALMANlearn2', 1, '(j)n', 'noise', [18 20 21 22]}, ... % # ,
    {'wh44wh39', 'RALMANlearn2', 1, '(d)kcc', 'noise', [18 20 21 22]}, ... % # 
    {'wh44wh39', 'RALMANlearn3', 1, '(m)d', 'noise', [17, 18, 20, 22]}, ... % #
    {'wh44wh39', 'RALMANlearn3', 1, '(n)h', 'noise', [18, 20, 22]}, ... % #  [18 20 22] ok
    {'wh44wh39', 'RALMANlearn3', 1, 'dk(c)c', 'noise'}, ... % # all
    {'wh44wh39', 'RALMANlearn3', 1, 'dkc(c)', 'noise', [17, 18, 20, 22]}, ... % # 
    {'wh44wh39', 'RALMANlearn3', 1, '(d)kcc', 'noise', [17, 18, 20, 22]}, ... % #
    {'wh44wh39', 'RALMANlearn3', 1, '(j)n', 'wn'}, ... % #
    {'wh44wh39', 'RALMANlearn3', 1, 'nh(h)', 'wn'}, ...
    {'wh44wh39', 'RALMANlearn4', 1, '(m)d', 'noise', [17 18 20 21 22]}, ... % @  ok
    {'wh44wh39', 'RALMANlearn4', 1, '(d)kc', 'noise', [17 18 20 21 22]}, ... % @ ok
    {'wh44wh39', 'RALMANlearn4', 1, 'dk(c)', 'noise'}, ... % @ r3emove
    {'wh44wh39', 'RALMANlearn4', 1, '(j)n', 'noise', [17 18 20 21 22]}, ... % @ ra fine
    {'wh44wh39', 'RALMANlearn4', 1, '(n)hh', 'noise', [17 18 20 21 22]}, ... % @# ra fine    {'wh44wh39', 'RALMANlearn4', 2, 'dkc(c)'}, ...
    {'wh72pk12', 'RALMANLearn3', 7, 'jrb(h)'}, ...
    {'wh72pk12', 'RALMANLearn3', 9, 'jrb(h)'}, ...
    {'wh72pk12', 'RALMANLearn3', 11, 'jrb(h)'}, ...
    {'wh72pk12', 'RALMANLearn4', 1, 'klb(h)', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn4', 1, '(g)', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn5', 1, 'iob(h)', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn5', 1, '(j)rb', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn5', 1, '(g)', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, 'klb(h)', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, '(g)', 'wn'}, ... % also noise.
    {'wh72pk12', 'RALMANLearn7', 1, '(j)rb', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, '(a)', 'wn'}, ... % also noise
    {'wh72pk12', 'RALMANLearn7', 1, '(j)kl', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, '(m)k', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, '(j)rb', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, 'j(r)b', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn7', 1, '(j)io', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn7', 2, 'klb(h)'}, ...
    {'wh72pk12', 'RALMANLearn7', 2, '(g)'}, ...
    {'wh72pk12', 'RALMANLearn7', 2, '(j)rb'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, '(m)k', 'wn'}, ... % noise too
    {'wh72pk12', 'RALMANLearn8', 1, 'm(k)', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, '(j)kl', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, 'j(k)l', 'wn'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, 'jrb(h)', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, '(j)rb', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, '(g)', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, 'kl(b)', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, 'klb(h)', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, '(j)io', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, 'io(b)', 'noise'}, ...
    {'wh72pk12', 'RALMANLearn8', 1, 'iob(h)', 'noise'}, ...
    {'gr48bu5', 'RALMANLearn2', 1, 'ab(h)', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn2', 1, '(a)j', 'wn'}, ... % @#
    {'gr48bu5', 'RALMANLearn2', 1, '(r)d', 'noise'}, ...
    {'gr48bu5', 'RALMANLearn2', 1, '(a)rd', 'noise'}, ...% @#
    {'gr48bu5', 'RALMANLearn2', 1, '(a)b', 'noise'}, ...% @#
    {'gr48bu5', 'RALMANLearn3', 1, 'jb(h)', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn3', 1, '(a)j', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn3', 1, '(a)rd', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn3', 1, '(a)bh', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn3', 1, '(r)rd', 'noise'}, ... % @# **********************
    {'gr48bu5', 'RALMANLearn3', 1, '(j)jb', 'noise'}, ... % @# ADDED NEW (3/17/19)
    {'gr48bu5', 'RALMANLearn4', 1, 'jb(h)', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn4', 1, '(a)rd', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn4', 1, '(a)bh', 'noise'}, ...
    {'gr48bu5', 'RALMANLearn4', 1, '(a)j', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn4', 1, 'a(b)h', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn4', 1, '(r)d', 'noise'}, ... % @# ADDED NEW (3/17)
    {'gr48bu5', 'RALMANLearn4', 1, '(r)rd', 'wn'}, ... % @# 
    {'gr48bu5', 'RALMANLearn5', 1, 'ab(h)', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn5', 1, '(a)j', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn5', 1, '(j)jb', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn5', 1, '(a)rd', 'wn'}, ... % @#
    {'gr48bu5', 'RALMANLearn5', 1, '(a)b', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn5', 1, '(r)rd', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn5', 1, '(r)d', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn5', 1, 'r(d)', 'noise'}, ... % @#
    {'gr48bu5', 'RALMANLearn6', 1, '(r)rd', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn6', 1, '(r)d', 'wn'}, ...
    {'gr48bu5', 'RALMANLearn6', 1, '(a)rd', 'wn'}, ... % @#
    {'gr48bu5', 'RALMANLearn6', 1, '(a)j', 'wn'}, ... % @#
    {'gr48bu5', 'RALMANLearn6', 1, '(a)bh', 'noise', [20 23]}, ... % @#
    {'gr48bu5', 'RALMANLearn6', 1, 'ab(h)', 'noise', [20 23]}, ...
    {'gr48bu5', 'RALMANLearn6', 1, 'jb(h)', 'wn'}, ...
    };


%     {'wh72pk12', 'RALMANLearn6', 1, '(j)io', 'wn'}, ... % ************************************

%% ========== fill in whether is "noise" or "wn"
% === if empty, then is "wn"
for i=1:length(SylsBad)
   if length(SylsBad{i})==4
       % then add on the end
      SylsBad{i} = [SylsBad{i} 'wn'];
   end
end


%% ========== ask whether the input matches any of the bad syls
sylbad = 0;
for j=1:length(SylsBad)
    tmp1 = strcmp(SylsBad{j}{1}, bname);
    tmp2 = strcmp(SylsBad{j}{2}, ename);
    tmp3 = SylsBad{j}{3}==swnum;
    tmp4 = strcmp(SylsBad{j}{4}, syltoken);
    
    % -- if this is chan to exclude, then abort
    if ~isempty(chanthis) & length(SylsBad{j})>5
        
        if any(ismember(SylsBad{j}{6}, chanthis))
        % then this channel should be excludede, i.e. this syl cannot be
        % bad for this chan
            continue
        end
    end
    
    if all([tmp1 tmp2 tmp3 tmp4]) & any(ismember(typesToRemove, SylsBad{j}{5}))
        sylbad=1;
        disp(['BAD SYL: ' bname '-' ename '-sw' num2str(swnum) '-' syltoken]);
        break
    end
end

