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


sylmax = max(DATBYREND.Sylcounter);


%% ==============
%% got thru all syls and ask whether they pass all criteria

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
%         snum = unique(DATBYREND.Sylnum(indsthis));
        
        % =============== passes syltype?        
        samethis = unique(DATBYREND.IsSame(indsthis));
        targthis = unique(DATBYREND.IsTarg(indsthis)); 
        if samethis~=getsame ...
                | targthis~=gettarg
            continue
        end
        
        % ============== bird to get?
        isSDP = unique(DATBYREND.Isfrom_SDP(indsthis));
        if ~isempty(birdstoget)
            
            if isnumeric(birdstoget)
                if ismember(bnum, birdstoget)==0
                    continue
                end
            elseif strcmp(birdstoget, 'notSDP')
                % then only if not from SDP experiments
                if isSDP==1
                    continue
                end
            elseif strcmp(birdstoget, 'onlySDP')
                if isSDP==0
                    continue
                end
            end
        end
        
        % ================ has baseline?
        numbaseinds = sum(DATBYREND.IsDurTrain(indsthis)==0);
        if numbaseinds < minbaserends
            continue
        end
        
        % ======= train rends
        numrends = sum(DATBYREND.IsDurTrain(indsthis)==1);
        if numrends < minrends
            continue
        end
        
        
        % === learning
        siglearn = unique(DATBYREND.SigLearn(indsthis));
        if getsiglearn==1
            if siglearn ==0
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


