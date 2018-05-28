function ffstruct = lt_neural_v2_CTXT_BRANCH_GetFF(analyfname, birdnum, neurnum,...
    branchnum, prms, classindsinorder)
%% optional
% classindsinorder = [2 4], then extracts just these classes. if ~exist,
% then extracts all classes, i.e. those that exist in the saved dir. NOTE:
% those classnums correspond to the classes in ALLBRANCH.FR, and those in
% CLASSES.

if ~exist('classindsinorder', 'var')
    classindsinorder = 1:1000; % basically all potentiall classes.
end

if isempty(classindsinorder)
    classindsinorder = 1:1000;
end

%%

% ============ load FF data
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

% ======= check that data exists
if exist([savedir '/' analyfname '/FF'], 'dir') ==0
    ffstruct = [];
    disp('no ffstruct found ..., retgurning ffstruct = []');
    return
end

% ============= filename is:
fname = [savedir '/' analyfname '/FF/' ...
    'bird' num2str(birdnum) '_neur' num2str(neurnum) '_branch' num2str(branchnum) ...
    '_classnum*.mat'];


%% iterate thru all classes and load
disp('---------');
classfnames = dir(fname);
count = 0;
for j=1:length(classfnames)
    
    % --------------------- make sure is one of desired classes
    indstr1 = strfind(classfnames(j).name, 'classnum');
    indstr2 = strfind(classfnames(j).name, '.mat');
    
    classthis = str2num(classfnames(j).name(indstr1+8:indstr2-1));
    if ~ismember(classthis, classindsinorder)
        continue
    end        
    
    % ---------------------- extract for this class.
    disp(classfnames(j).name);
    tmp = load([savedir '/' analyfname '/FF/' classfnames(j).name]);
    t_ff = tmp.t_ff; % [t, ff]
    
    % --- ignore if smaller than sample
    if size(t_ff,1)<prms.minN
        continue
    end
    
    % ==== output
    count = count+1;
    ffstruct.classnum(count).t_ff = t_ff;
    ffstruct.classnum(count).name = classfnames(j).name;
    ffstruct.classnum(count).ind_classorig = classthis;
end

%% ======== check whether match regexp
if (0)
    if birdnum>1
        tmp = load('AllBirdNeurBranchClass.mat');
        indtmp = tmp.AllBirdNeurBranchClass(:, 1)==birdnum & ...
            tmp.AllBirdNeurBranchClass(:,2)==neurnum ...
            & tmp.AllBirdNeurBranchClass(:,3)==branchnum;
        
        assert(sum(indtmp) == length(ffstruct.classnum), 'sadfasdf');
    end
end
