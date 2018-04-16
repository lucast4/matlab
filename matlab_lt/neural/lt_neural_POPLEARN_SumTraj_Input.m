%% input data, hand entered
% ========== each linked by hand to specific MOTIFSTATS_Compiled structure

%% === only pu69
clear metadatcell;
clear metadatstruct;
% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};

ind = 1;
metadatstruct(ind).bname = 'pu69wh78';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'aa(b)';
metadatstruct(ind).base_set = [1 2];
metadatstruct(ind).base_switches = {'', '1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [2];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{}};
metadatstruct(ind).targlearndir = 1;


%% only good learning expts
if (0)
   load  MOTIFSTATS_Compiled_31Mar2018_1639;
end
clear metadatcell;
clear metadatstruct;
% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};
ind = 1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'nh(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [4 5 6 7 8];
metadatstruct(ind).train_switches = {'', '', '', '', ''};
metadatstruct(ind).train_sylsignore = {{}, {}, {}, {}, {}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'n(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [1 2 3];
metadatstruct(ind).train_switches = {'', '', '1r'};
metadatstruct(ind).train_sylsignore = {{'nh(h)'}, {'nh(h)'}, {'nh(h)'}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn2';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [1];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{'cb(b)'}};
metadatstruct(ind).train_set = [1];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [3 2 1];
metadatstruct(ind).train_switches = {'', '7l', '9l'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}, {'cb(b)'}, {'cb(b)'}};
metadatstruct(ind).targlearndir = -1;


ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'cb(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [5 4 2];
metadatstruct(ind).train_switches = {'', '', '7r'};
metadatstruct(ind).train_sylsignore = {{}, {}, {'c(b)'}};
metadatstruct(ind).targlearndir = -1;


% ======================== PU69
ind = ind+1;
metadatstruct(ind).bname = 'pu69wh78';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'aa(b)';
metadatstruct(ind).base_set = [1 2];
metadatstruct(ind).base_switches = {'', '1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [2];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{}};
metadatstruct(ind).targlearndir = 1;

% ind = ind+1;
% metadatstruct(ind).bname = 'pu69wh78';
% metadatstruct(ind).exptname = 'RALMANlearn2';
% metadatstruct(ind).targsyl = 'aa(b)';
% metadatstruct(ind).base_set = [1];
% metadatstruct(ind).base_switches = {'1l'};
% metadatstruct(ind).base_sylsignore = {{}};
% metadatstruct(ind).train_set = [1];
% metadatstruct(ind).train_switches = {'1r'};
% metadatstruct(ind).train_sylsignore = {{}};
% metadatstruct(ind).targlearndir = 1;

% ind = ind+1;
% metadatstruct(ind).bname = 'pu69wh78';
% metadatstruct(ind).exptname = 'RALMANOvernightLearn1';
% metadatstruct(ind).targsyl = 'aab(h)';
% metadatstruct(ind).base_set = [4];
% metadatstruct(ind).base_switches = {''};
% metadatstruct(ind).base_sylsignore = {{}};
% metadatstruct(ind).train_set = [3];
% metadatstruct(ind).train_switches = {''};
% metadatstruct(ind).train_sylsignore = {{}};
% metadatstruct(ind).targlearndir = 1;

%% ========== wh44/pu69 all, including all epochs.
% regardless of whether good learning. except cb(b) in first epoch for 
% RALMANlearn1, since no learning at all.
if (0)
   load  MOTIFSTATS_Compiled_31Mar2018_1639;
end
clear metadatcell;
clear metadatstruct;
% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};
ind = 1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'nh(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [4 5 6 7 8];
metadatstruct(ind).train_switches = {'', '', '', '', ''};
metadatstruct(ind).train_sylsignore = {{}, {}, {}, {}, {}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'n(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [1 2 3];
metadatstruct(ind).train_switches = {'', '', '1r'};
metadatstruct(ind).train_sylsignore = {{'nh(h)'}, {'nh(h)'}, {'nh(h)'}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn2';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [1];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{'cb(b)'}};
metadatstruct(ind).train_set = [1];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [6 3 2 1];
metadatstruct(ind).train_switches = {'2r', '', '7l', '9l'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}, {'cb(b)'}, {'cb(b)'}, {'cb(b)'}};
metadatstruct(ind).targlearndir = -1;


ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'cb(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [5 4 2];
metadatstruct(ind).train_switches = {'', '', '7r'};
metadatstruct(ind).train_sylsignore = {{}, {}, {'c(b)'}};
metadatstruct(ind).targlearndir = -1;


% ======================== PU69
ind = ind+1;
metadatstruct(ind).bname = 'pu69wh78';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'aa(b)';
metadatstruct(ind).base_set = [1 2];
metadatstruct(ind).base_switches = {'', '1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [2];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'pu69wh78';
metadatstruct(ind).exptname = 'RALMANlearn2';
metadatstruct(ind).targsyl = 'aa(b)';
metadatstruct(ind).base_set = [1];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [1];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{}};
metadatstruct(ind).targlearndir = 1;

ind = ind+1;
metadatstruct(ind).bname = 'pu69wh78';
metadatstruct(ind).exptname = 'RALMANOvernightLearn1';
metadatstruct(ind).targsyl = 'aab(h)';
metadatstruct(ind).base_set = [4];
metadatstruct(ind).base_switches = {''};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [3];
metadatstruct(ind).train_switches = {''};
metadatstruct(ind).train_sylsignore = {{}};
metadatstruct(ind).targlearndir = 1;

%% ========== 3/31/18 - wh44 all, including all epochs.
% regardless of whether good learning. except cb(b) in first epoch for 
% RALMANlearn1, since no learning at all.
if (0)
    load MOTIFSTATS_Compiled_30Mar2018_1422
end

clear metadatcell;
clear metadatstruct;
% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};
ind = 1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'nh(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [4 5 6 7 8];
metadatstruct(ind).train_switches = {'', '', '', '', ''};
metadatstruct(ind).train_sylsignore = {{}, {}, {}, {}, {}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'n(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [1 2 3];
metadatstruct(ind).train_switches = {'', '', '1r'};
metadatstruct(ind).train_sylsignore = {{'nh(h)'}, {'nh(h)'}, {'nh(h)'}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn2';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [1];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{'cb(b)'}};
metadatstruct(ind).train_set = [1];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [6 3 2 1];
metadatstruct(ind).train_switches = {'2r', '', '7l', '9l'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}, {'cb(b)'}, {'cb(b)'}, {'cb(b)'}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'cb(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [5 4 2];
metadatstruct(ind).train_switches = {'', '', '7r'};
metadatstruct(ind).train_sylsignore = {{}, {}, {'c(b)'}};



%% ========== 3/31/18 - only first day

clear metadatcell;
clear metadatstruct;
% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};
% ind = 1;
% metadatstruct(ind).bname = 'wh44wh39';
% metadatstruct(ind).exptname = 'RALMANlearn3';
% metadatstruct(ind).targsyl = 'nh(h)';
% metadatstruct(ind).base_set = [3];
% metadatstruct(ind).base_switches = {'1l'};
% metadatstruct(ind).base_sylsignore = {{}};
% metadatstruct(ind).train_set = [4 5 6 7 8];
% metadatstruct(ind).train_switches = {'', '', '', '', ''};
% metadatstruct(ind).train_sylsignore = {{}, {}, {}, {}, {}};

ind = 1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn3';
metadatstruct(ind).targsyl = 'n(h)';
metadatstruct(ind).base_set = [3];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [1 2 3];
metadatstruct(ind).train_switches = {'', '', '1r'};
metadatstruct(ind).train_sylsignore = {{'nh(h)'}, {'nh(h)'}, {'nh(h)'}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn2';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [1];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{'cb(b)'}};
metadatstruct(ind).train_set = [1];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'c(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [6];
metadatstruct(ind).train_switches = {'2r'};
metadatstruct(ind).train_sylsignore = {{'cb(b)'}};

ind = ind+1;
metadatstruct(ind).bname = 'wh44wh39';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'cb(b)';
metadatstruct(ind).base_set = [6];
metadatstruct(ind).base_switches = {'1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [6];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{}};


%% ===== just pu69
if (0)
    load MOTIFSTATS_Compiled_22Feb2018_1202;
end

clear metadatcell;
clear metadatstruct;
% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};
% ind = 1;
% metadatstruct(ind).bname = 'wh44wh39';
% metadatstruct(ind).exptname = 'RALMANlearn3';
% metadatstruct(ind).targsyl = 'nh(h)';
% metadatstruct(ind).base_set = [3];
% metadatstruct(ind).base_switches = {'1l'};
% metadatstruct(ind).base_sylsignore = {{}};
% metadatstruct(ind).train_set = [4 5 6 7 8];
% metadatstruct(ind).train_switches = {'', '', '', '', ''};
% metadatstruct(ind).train_sylsignore = {{}, {}, {}, {}, {}};

ind = 1;
metadatstruct(ind).bname = 'pu69wh78';
metadatstruct(ind).exptname = 'RALMANlearn1';
metadatstruct(ind).targsyl = 'aa(b)';
metadatstruct(ind).base_set = [1 2];
metadatstruct(ind).base_switches = {'', '1l'};
metadatstruct(ind).base_sylsignore = {{}};
metadatstruct(ind).train_set = [2];
metadatstruct(ind).train_switches = {'1r'};
metadatstruct(ind).train_sylsignore = {{'nh(h)'}};


%% old version, using cell array
% metadatcell{1} = {'wh44wh39', 'RALMANlearn3', 'nh(h)', [3], {{}}, ...
%     [4 5 6 7 8], {{}, {}, {}, {}, {}}};
% metadatcell{2} = {'wh44wh39', 'RALMANlearn3', 'n(h)', [3], {{}}, ...
%     [1 2], {{'nh(h)'}, {'nh(h)'}}};
% metadatcell{3} = {'wh44wh39', 'RALMANlearn2', 'c(b)', [1], {{'cb(b)'}}, ...
%     [1], {{'cb(b)'}}};
% metadatcell{4} = {'wh44wh39', 'RALMANlearn1', 'c(b)', [6], {{}}, ...
%     [6 3 2 1], {{'cb(b)'}, {'cb(b)'}, {'cb(b)'}, {'cb(b)'}}};
% metadatcell{5} = {'wh44wh39', 'RALMANlearn1', 'cb(b)', [6], {{}}, ...
%     [6 5 4 2], {{}, {}, {}, {'c(b)'}}};
