function [indmotif_all, motiflist_out] = lt_neural_QUICK_MotifID(birdname, motifregexp)
%% lt 10/14/18 - tells you position of a given motif name in global rference frame
% NOTE: if multiple motifs spelled differently, but refer to same thing
% (e.g. a(b)h and a(b)hh), then can make them the same below.

%  input
% motifregexp = 'a(b)'; CAN ALSO BE cell array Nx1. e.g.: MUST ALL BE SAME
% BIRD ...
% motifs =
%
%   6×1 cell array
%
%     '(j)jbhh'
%     '(j)jbhh'
%     '(j)jbhh'
%     '(j)jbhh'
%     'j(j)bhh'
%     'j(j)bhh'
%
% birdname = 'wh44wh39';

% indmotif is ind within entire lsit of motifs.

% motiflist_out, for this bird, e.g. for labeling stuff.

%%

if ~exist('motifregexp', 'var')
    motifregexp = '';
end

%% TOOL TO CHECK WHAT MOTIFS EXIST FOR A GIVEN DATA STRUCTURE
if (0)
    % ==== 1) quickly list all motifs
    disp('========================');
    numbirds = length(MOTIFSTATS_Compiled.birds);
    for i=1:numbirds
        bname = MOTIFSTATS_Compiled.birds(i).birdname;
        MotiflistAll = {};
        numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
        for ii=1:numneurons
            
            motiflist = [MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str];
            MotiflistAll = [MotiflistAll motiflist];
            
        end
        MotiflistAll = unique(MotiflistAll);
        disp(['BIRD: ' birdname]);
        disp(MotiflistAll);
    end
    
end

%% ========================= PARAMS, hand enter these. they are in order of singing.
MotifDatabase = struct;

% =====================================
MotifDatabase.bird(1).birdname = 'pu69wh78';
MotifDatabase.bird(1).motifs = {...
    {'(a)ab'}, ...
    {'a(a)b'}, ...
    {'aa(b)'}, ...
    {'aab(h)'}, ...
    {'aabh(h)'}, ...
    {'h(g)'}, ...
    {'(j)jb', '(j)jbhh'}, ...
    {'j(j)b', 'j(j)bhh'}, ...
    {'jj(b)', 'jj(b)hh'}, ...
    {'jjb(h)', 'jjb(h)h'}, ...
    {'jjbh(h)'}, ...
    {'jjbhh(h)'}, ...
    {'jjbhh(g)'}, ...
    };

MotifDatabase.bird(2).birdname = 'wh44wh39';
MotifDatabase.bird(2).motifs = {...
    { '(m)d' }, ...
    {'(d)kc', '(d)kcc'}, ...
    {'d(k)c', 'd(k)cc'}, ...
    {'dk(c)', 'dk(c)c'}, ...
    {'dkc(c)'}, ...
    {'c(b)'}, ...
    {'cb(b)'}, ...
    {'(j)n'}, ...
    {'(n)h', '(n)hh'}, ...
    {'n(h)', 'n(h)h'}, ...
    {'nh(h)'}, ...
    };


%% ============= CONVERT EACH TO MOTIF LISTS (removing redundant - only for output)

for i=1:length(MotifDatabase.bird)
    
    motiflist_short = {};
    for ii=1:length(MotifDatabase.bird(i).motifs)
       
        motifthis = MotifDatabase.bird(i).motifs{ii}{1};
        motiflist_short = [motiflist_short motifthis];
        
        
    end
    
    MotifDatabase.bird(i).motiflist_short = motiflist_short;
end

%%

indbird = strcmp({MotifDatabase.bird.birdname}, birdname);
assert(sum(indbird)==1, 'PROBLEM - this bird not entered yet');

motiflist_out = MotifDatabase.bird(indbird).motiflist_short;

if isempty(motifregexp)
    indmotif_all = [];
else
    if iscell(motifregexp)
    indmotif_all = [];
    for j=1:length(motifregexp)
        motifthis = motifregexp{j};
        functmp = @(X)any(strcmp(X, motifthis));
        indmotif = find(cellfun(functmp, MotifDatabase.bird(indbird).motifs));
        assert(length(indmotif)<2, 'PROBLEM - entered same motif in multiple entries ...');
        if isempty(indmotif)
            keyboard
            disp('PROBLEM --- this motif not entered for this bird...');
        end
        indmotif_all = [indmotif_all; indmotif];
    end
else % then is a single motif
    functmp = @(X)any(strcmp(X, motifregexp));
    indmotif_all = find(cellfun(functmp, MotifDatabase.bird(indbird).motifs));
    assert(length(indmotif_all)<2, 'PROBLEM - entered same motif in multiple entries ...');
    if isempty(indmotif_all)
        keyboard
        disp('PROBLEM --- this motif not entered for this bird...');
    end
    end
end
