function OUTSTRUCT_XCOV = lt_neural_POPLEARN_XCovExtr_post(OUTSTRUCT_XCOV, ...
    SwitchCohStruct, SummaryStruct, SwitchStruct)

%%


% ======== name change, so OUTSTRUCT_XCOV matches OUTSTRUCT
OUTSTRUCT_XCOV.switch = OUTSTRUCT_XCOV.swnum;
OUTSTRUCT_XCOV.chanpair = OUTSTRUCT_XCOV.neurpair;



% ####################### EXTRACT CERTAIN THINGS INTO OUTSTRUCT_XCOV
% ====== 1) MOTIFNAME
motifname_all = {};
chanpair_actual = [];
for i=1:length(OUTSTRUCT_XCOV.bnum)
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    swnum = OUTSTRUCT_XCOV.swnum(i);
    motifnum = OUTSTRUCT_XCOV.motifnum(i);

    motifname = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swnum).motifnum(motifnum).motifname;
    motifname_all = [motifname_all; motifname];
    
    chanpair_actual = [chanpair_actual; ...
        [SummaryStruct.birds(bnum).neurons(OUTSTRUCT_XCOV.neurpair(i,:)).channel]];
end

OUTSTRUCT_XCOV.motifname = motifname_all;
OUTSTRUCT_XCOV.chanpair_actual = chanpair_actual;

% ======



% ======================= [EXTRACT] unique motif ID
allmotifID = [];
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    bnum = OUTSTRUCT_XCOV.bnum(i);
    motifname = OUTSTRUCT_XCOV.motifname{i};
    
    bname = SwitchStruct.bird(bnum).birdname;
    
    motifID = lt_neural_QUICK_MotifID(bname, motifname);
    
    allmotifID = [allmotifID; motifID];    
    
end
OUTSTRUCT_XCOV.motifID_unique = allmotifID;