function [XDAT] = lt_neural_Coher_XDATSET_Extract(savedir, indtoget_b_e_s, ...
    SwitchStruct)

%% lt 1/3/18 - extracts previously saveld across-data-set stuff


% ==== COLLECT LIST OF ALL BIRDS, EXPT, and SW
fnames = dir([savedir '/*.mat']);

all_bname = {};
all_ename = {};
all_swnum = [];
all_data = [];
all_savemarker = {};
for i=1:length(fnames)
    fthis = fnames(i).name;
    
    % --- from file name, dewtermint things
    indtmp = strfind(fthis, '_');
    bname = fthis(1:indtmp(1)-1);
    ename = fthis(indtmp(1)+1:indtmp(2)-1);
    swnum = str2num(fthis(indtmp(2)+3:indtmp(3)-1));
    savemarker = fthis(indtmp(3)+1:indtmp(3)+14);
    %    swnum = fthis(indtmp(2)+3:indtmp(3)-1)
    
    % ========== DECIDE IF SKIP THIS CASE
    if ~isempty(indtoget_b_e_s)
        bnum = find(strcmp({SwitchStruct.bird.birdname}, bname));
        enum = find(strcmp({SwitchStruct.bird(bnum).exptnum.exptname}, ename));
        
        % ==== decide whether to skip
            if ~ismember([bnum enum swnum], indtoget_b_e_s, 'rows')
                continue
            end
    end
    
    all_bname = [all_bname; bname];
    all_ename = [all_ename; ename];
    all_swnum = [all_swnum; swnum];
    all_savemarker = [all_savemarker; savemarker];
    
    % --- extract data
    tmp = load([savedir '/' fthis]);
    all_data = [all_data; tmp.savestruct];
    
end

%% ==== output
XDAT.all_bname = all_bname;
XDAT.all_ename = all_ename;
XDAT.all_swnum = all_swnum;
XDAT.all_savemarker = all_savemarker;
XDAT.all_data = all_data;
