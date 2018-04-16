%% lt 4/4/18 - looking at timecourse of generalization
%% ################## COLLECT ALL LEARNING DATA
GenStruct = struct;
ind = 0;

% ================================= SEQ DEP PITCH EXPERIMENTS



% ================================= NEURAL RECORDING EXPERIMENTS


% ======= pu69wh78 - RALMANlearn1
birdname = 'pu69wh78';
exptname = 'RALMANlearn1';
[DATSTRUCT, TargSyl, MotifList] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).targsyl = TargSyl;
GenStruct.expt(ind).motiflist= MotifList;



%% STORE HERE THE PARAMS USED IN EACH EXPERIMENT





%% PLOT EACH EXPERIMENT (OVERLAY TARGET VS. NONTARGETS)



%% 


