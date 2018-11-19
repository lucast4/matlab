function OUTSTRUCT = lt_neural_LFP_GetLearnDir(OUTSTRUCT, SwitchStruct)
%%  lt 11/11/18 - gets 1) learn dir, and 2) ff diff

OUTSTRUCT.learndirTarg = [];
OUTSTRUCT.ff_WNminusBase = [];
for i=1:length(OUTSTRUCT.bnum)
    
    bnum = OUTSTRUCT.bnum(i);
    enum = OUTSTRUCT.enum(i);
    sw = OUTSTRUCT.switch(i);
    %     mm = OUTSTRUCT.motifnum(i);
    
    ldirs = cell2mat(SwitchStruct.bird(bnum).exptnum(enum).switchlist(sw).learningDirs(2:2:end));
    learndir = unique(ldirs);
    
    OUTSTRUCT.learndirTarg = [OUTSTRUCT.learndirTarg; learndir];
    
    % ---- ff diff
    ff = OUTSTRUCT.ffvals{i};
    inds_base = OUTSTRUCT.indsbase_epoch{i};
    inds_WN = OUTSTRUCT.indsWN_epoch{i};
    
    ffdiff = mean(ff(inds_WN)) - mean(ff(inds_base));
    OUTSTRUCT.ff_WNminusBase = [OUTSTRUCT.ff_WNminusBase; ffdiff];
end
