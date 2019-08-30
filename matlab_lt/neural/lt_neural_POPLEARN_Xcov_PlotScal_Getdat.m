function lt_neural_POPLEARN_Xcov_PlotScal_Getdat()
%% lt 6/2/19 - gets all relevant data into cell arrays
% all datapoints (even low level: syl/chan)


Yall_scal = {};
Yall_bnum = [];
Yall_enum = [];
Yall_chanID = [];
Yall_istarg = {};
Yall_issame = {};
Yall_sw = [];
for i=1:length(indsgrp_chanpair_unique)
    
    indsthis = indsgrp_chanpair == indsgrp_chanpair_unique(i);
    
    bnum = unique(OUTSTRUCT.bnum(indsthis));
    enum = unique(OUTSTRUCT.enum(indsthis));
    swnum = unique(OUTSTRUCT.switch(indsthis));
    
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    % =========================== DATA
    cohscal = Yscalar(indsthis);
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    chnums = OUTSTRUCT.chanpair_actual(indsthis,:);
    neurpair = OUTSTRUCT.neurpair(indsthis,:);
    
    % ================================ SAVE OUTPUT
    Yall_scal = [Yall_scal; cohscal];
    Yall_bnum = [Yall_bnum ; bnum];
    Yall_enum = [Yall_enum ; enum];
    Yall_chanID = [Yall_chanID; indsgrp_chanpair_unique(i)];
    Yall_istarg = [Yall_istarg; istarg];
    Yall_issame = [Yall_issame; issame];
    Yall_sw = [Yall_sw; swnum];
end
