function [Inds_sylcounter, Inds_birdnum, Inds_exptnum, ...
    IsSame, IsTarg] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
    syltype, getsiglearn, birdstoget, minbaserends, minrends)

%% 8/2/18 - lt, pull out inds of syls to plot that pass criteria entered.


%% params
% minrendsinbin = 4;
%
% minbaserends = 100; % has at least this many baseline renditions.
% getsiglearn = 0;
%
% % birdstoget = 13:17;
% % birdstoget = 1;
% % birdstoget = 'notSDP'; % e.g. generalization struct
% % birdstoget = 'onlySDP'; % only from SDP experiments
% % birdstoget = 'all_preferNoSDP' % all, but if multipel version of same
% expt, then takes noSDP
% leave empty to get all birds

% % syltype = 'same';
% syltype = 'same';
% % syltype = 'diff';
%


%%

if strcmp(syltype, 'same')
    getsame = 1;
    gettarg = 0;
elseif strcmp(syltype, 'diff')
    getsame = 0;
    gettarg = 0;
elseif strcmp(syltype, 'targ')
    getsame = 1;
    gettarg = 1;
elseif strcmp(syltype, 'all')
    getsame = [];
    gettarg = [];
elseif strcmp(syltype, 'nontarg')
    getsame = [];
    gettarg = 0;
end




%% if want to remove redundant cases, then do that here
if strcmp(birdstoget, 'all_preferNoSDP') | strcmp(birdstoget, 'all_preferSDP')
    
    if length(DATBYREND.IsWN) ~= length(DATBYREND.Birdnum)
        % not sure why, just remove this field
        DATBYREND = rmfield(DATBYREND, 'IsWN');
    end
    
    nbirds = length(TrialStruct.birds);
    for i=1:nbirds
        numexpt = length(TrialStruct.birds(i).exptnum);
        exptbad = [];
        
        % --- go thru all pairs of expts. flag if idnetical
        for j=1:numexpt
            for jj=j+1:numexpt
                ename1 = TrialStruct.birds(i).exptnum(j).exptname;
                ename2 = TrialStruct.birds(i).exptnum(jj).exptname;
                
                if strcmp(ename1, ename2)
                    % ---- figure out which one is SDP
                    sdp1 = logical(unique(DATBYREND.Isfrom_SDP(DATBYREND.Birdnum==i & DATBYREND.Exptnum==j)));
                    sdp2 = logical(unique(DATBYREND.Isfrom_SDP(DATBYREND.Birdnum==i & DATBYREND.Exptnum==jj)));
                    
                    if isempty(sdp1) | isempty(sdp2)
                        % --- then fine, only one has any extracted data
                        continue
                    else
                        assert(~(sdp1 & sdp2), 'problem, why both are sdp?')
                        assert(sdp1 | sdp2, 'problem, why neither are sdp?')
                        
                        % --- not fine - remove the one that you don't want
                        if strcmp(birdstoget, 'all_preferNoSDP')
                            % -- then remove SDP
                            if sdp1==1
                                indbad = j;
                            else
                                indbad = jj;
                            end
                        elseif strcmp(birdstoget, 'all_preferSDP')
                            if sdp1==1
                                indbad = jj;
                            else
                                indbad = j;
                            end
                        end
                    end
                    
                    exptbad = [exptbad indbad];
                end
            end
        end
        
        % ======== DO removal
        for j=exptbad
            indsgood = ~(DATBYREND.Birdnum==i & DATBYREND.Exptnum==j);
            DATBYREND = lt_structure_subsample_all_fields(DATBYREND, indsgood, 1);
        end
    end
end
%% ==============
%% got thru all syls and ask whether they pass all criteria
sylmax = max(DATBYREND.Sylcounter);

Inds_sylcounter = [];
Inds_birdnum = [];
Inds_exptnum = [];
IsSame = [];
IsTarg = [];

for ss = 1:sylmax
    
    indsthis = DATBYREND.Sylcounter==ss;
    
    if ~any(indsthis)
        continue
    end
    
    % ---- want this syl?
    bnum = unique(DATBYREND.Birdnum(indsthis));
    enum = unique(DATBYREND.Exptnum(indsthis));
    snum = unique(DATBYREND.Sylnum(indsthis));
    
    % =============== passes syltype?
    samethis = unique(DATBYREND.IsSame(indsthis));
    targthis = unique(DATBYREND.IsTarg(indsthis));
    if samethis~=getsame ...
            | targthis~=gettarg
        disp('removing, syl type bad');
        continue
    end
    
    % ============== bird to get?
    isSDP = unique(DATBYREND.Isfrom_SDP(indsthis));
    if ~isempty(birdstoget)
        
        if isnumeric(birdstoget)
            if ismember(bnum, birdstoget)==0
                disp('removing');
                continue
            end
        elseif strcmp(birdstoget, 'notSDP')
            % then only if not from SDP experiments
            if isSDP==1
                disp('removing');
                continue
            end
        elseif strcmp(birdstoget, 'onlySDP')
            if isSDP==0
                disp('removing');
                continue
            end
        elseif strcmp(birdstoget, 'all_preferNoSDP') | strcmp(birdstoget, 'all_preferSDP')
            % ignore, already took care of remove redundancies above.
        end
    end
    
    % ================ has baseline?
    numbaseinds = sum(DATBYREND.IsDurTrain(indsthis)==0);
    if numbaseinds < minbaserends
        disp('removing, baseline rends not enough');
        continue
    end
    
    % ======= train rends
    numrends = sum(DATBYREND.IsDurTrain(indsthis)==1);
    if numrends < minrends
        bname = TrialStruct.birds(bnum).birdname;
        ename = TrialStruct.birds(bnum).exptnum(enum).exptname;
        sylname = TrialStruct.birds(bnum).exptnum(enum).sylnum(snum).syl;
        disp(['REMOVING becuase low rends: ' bname '-' ename '-' sylname]);
        continue
    end
    
    
    % === learning
    siglearn = unique(DATBYREND.SigLearn(indsthis));
    if getsiglearn==1
        if siglearn ==0
            disp('removing, poor learning');
            continue
        end
    end
    
    % ================ COLLECT IF UP TO HERE AND HAVE NOT CONTINUED
    Inds_sylcounter = [Inds_sylcounter ss];
    Inds_birdnum = [Inds_birdnum bnum];
    Inds_exptnum = [Inds_exptnum enum];
    IsSame = [IsSame samethis];
    IsTarg = [IsTarg targthis];
end


