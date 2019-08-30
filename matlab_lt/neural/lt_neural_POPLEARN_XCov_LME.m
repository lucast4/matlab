function lt_neural_POPLEARN_XCov_LME(Yscalar, OUTSTRUCT_XCOV, ...
    SwitchStruct, MOTIFSTATS_Compiled, PARAMS, clim, ...
    onlygoodexpt, expttype, onlyIfSameType)

%%
disp('ONLY LOOKING AT FIRST SWITCH (need to modify code if want otherwise)!');
disp('ASSUMES THAT ALL data for a given b/e are from one switch.. [must change if not true]');

pause;

%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter outstruct
    [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(...
        OUTSTRUCT_XCOV, SwitchStruct, expttype);
    % ===== filter the scalar array
    Yscalar = Yscalar(indstokeep);
end

% ====== SANITY CHECK.
assert(length(Yscalar) == length(OUTSTRUCT_XCOV.bnum));


%% ===== only same type?
if onlyIfSameType==1
    [indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch});
    indstokeep = [];
    
    for j=1:length(indsgrp_switch_unique)
        
        indsthis = find(indsgrp_switch==j);
        
        if any(OUTSTRUCT_XCOV.issame(indsthis)==1)
            indstokeep = [indstokeep; indsthis];
        end
    end
    
    OUTSTRUCT_XCOV = lt_structure_subsample_all_fields(OUTSTRUCT_XCOV, indstokeep, 1);
    Yscalar = Yscalar(indstokeep);
end

%% ============= syllable type, into categories

[indssyltype, ~, xcell] = lt_tools_grp2idx({OUTSTRUCT_XCOV.issame, OUTSTRUCT_XCOV.istarg});

% 10 3 (same)
% 00 1 (diff)
% 01 2 (targ)

%% ============= 1) GET UNIQUE INDS FOR EACH SWITCH, AND EACH CHAN PAIR
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch, ...
    OUTSTRUCT_XCOV.chanpair});

indsexpt = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum});



% =============== [EXTRACT DATA]
dat = table(Yscalar, categorical(OUTSTRUCT_XCOV.bnum), categorical(indsexpt), ...
    categorical(indsgrp_chanpair), OUTSTRUCT_XCOV.istarg, OUTSTRUCT_XCOV.issame,...
    categorical(indssyltype), 'VariableNames', {'Yresponse', 'bnum', 'exptnum', 'chanpairID', 'istarg', 'issame', 'syltype'});




%% ============== model [PLAYING AROUND - REAL ONE BELOW]

% formula = 'Yresponse ~ istarg + (1|chanpairID) + (-1 + istarg|chanpairID)';
% formula = 'Yresponse ~ istarg + (istarg|chanpairID)';
% formula = 'Yresponse ~ istarg + (istarg|bnum) + (istarg|bnum:exptnum)';
% formula = 'Yresponse ~ istarg + (istarg|bnum) + (istarg|bnum:exptnum) + (istarg|bnum:exptnum:chanpairID)';
formula = 'Yresponse ~ istarg + (istarg|exptnum) + (istarg|exptnum:chanpairID)';
% formula = 'Yresponse ~ istarg + (1|exptnum) + (-1+exptnum|exptnum) + (istarg|exptnum:chanpairID)';

lme = fitlme(dat, formula, 'StartMethod', 'random');

formula = 'Yresponse ~ 1 + (istarg|bnum) + (istarg|bnum:exptnum) + (istarg|bnum:exptnum:chanpairID)';
lmesimple = fitlme(dat, formula, 'StartMethod', 'random');

%% =========== simulated likelihood test.
% [results, siminfo] = compare(lmesimple, lme, 'CheckNesting', 1, 'Nsim', 2);


%% ============== model

% formula = 'Yresponse ~ istarg + (1|chanpairID) + (-1 + istarg|chanpairID)';

if onlyIfSameType==1
    % --- use categories for the three syl types
    formula = 'Yresponse ~ syltype + (istarg|exptnum) + (istarg|exptnum:chanpairID)';
else
    formula = 'Yresponse ~ istarg + (istarg|exptnum) + (istarg|exptnum:chanpairID)';
end
% formula = 'Yresponse ~ istarg + (istarg|bnum) + (istarg|bnum:exptnum) + (istarg|bnum:exptnum:chanpairID)';

lme = fitlme(dat, formula, 'StartMethod', 'random');


% ################## plot
lt_neural_POPLEARN_XCOV_LME_subplot;


%% ============== model [just targ syls, change from baseline?

% formula = 'Yresponse ~ istarg + (1|chanpairID) + (-1 + istarg|chanpairID)';


formula = 'Yresponse ~ 1 + (1|exptnum)';
% formula = 'Yresponse ~ istarg + (istarg|bnum) + (istarg|bnum:exptnum) + (istarg|bnum:exptnum:chanpairID)';

dat = dat(dat.istarg==1,:);
lme = fitlme(dat, formula, 'StartMethod', 'random');


% ################## plot
lt_neural_POPLEARN_XCOV_LME_subplot;

% ============== COEFFICIENT TEST?
% coefTest(lme, [0 1]);
% anova(lme, 'DFMethod', 'none')

