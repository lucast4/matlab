function [indmotif_all, motiflist_out, motifposition_all] = ...
    lt_neural_QUICK_MotifID(birdname, motifregexp)
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

% motifposition = [1 1; 1 2], i.e. [motif number, and position within
% motif.

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
        disp(['BIRD: ' bname]);
        disp(MotiflistAll);
    end
end


%% ========================= PARAMS, hand enter these. they are in order of singing.
MotifDatabase = struct;

% =====================================
MotifDatabase.bird(1).birdname = 'pu69wh78';
MotifDatabase.bird(1).motifs = {...
    {'(a)ab', 1}, ...
    {'a(a)b', 1}, ...
    {'aa(b)', 1}, ...
    {'aab(h)', 1}, ...
    {'aabh(h)', 1}, ...
    {'h(g)', 1}, ...
    {'(j)jb', '(j)jbhh', 2}, ...
    {'j(j)b', 'j(j)bhh', 2}, ...
    {'jj(b)', 'jj(b)hh', 2}, ...
    {'jjb(h)', 'jjb(h)h', 2}, ...
    {'jjbh(h)', 2}, ...
    {'jjbhh(h)', 2}, ...
    {'jjbhh(g)', 2}, ...
    };

MotifDatabase.bird(2).birdname = 'wh44wh39';
MotifDatabase.bird(2).motifs = {...
    { '(m)d' , 1}, ...
    {'(d)kc', '(d)kcc', 1}, ...
    {'d(k)c', 'd(k)cc', 1}, ...
    {'dk(c)', 'dk(c)c', 1}, ...
    {'dkc(c)', 1}, ...
    {'c(b)', 1}, ...
    {'cb(b)', 1}, ...
    {'(j)n', 2}, ...
    {'(n)h', '(n)hh', 2}, ...
    {'n(h)', 'n(h)h', 2}, ...
    {'nh(h)', 2}, ...
    };

MotifDatabase.bird(3).birdname = 'wh72pk12';
MotifDatabase.bird(3).motifs = {...
    {'(j)rb', 1}, ...
    {'j(r)b', 1}, ...
    {'jr(b)', 1}, ...
    {'jrb(h)', 1}, ...
    {'(j)kl', 2}, ...
    {'j(k)l', 2}, ...
    {'jk(l)', 2}, ...
    {'kl(b)', 2}, ...
    {'klb(h)', 2}, ...
    {'(m)k', 3}, ...
    {'m(k)', 3}, ...
    {'mk(d)', 3}, ...
    {'(j)io', 4}, ...
    {'(i)ob', 4}, ...
    {'i(o)b', 4}, ...
    {'io(b)', 4}, ...
    {'iob(h)', 4}, ...
    {'(a)', 5}, ...
    {'(g)', 5}};


MotifDatabase.bird(4).birdname = 'gr48bu5';
MotifDatabase.bird(4).motifs = {...
    {'(r)rd', 1}, ...
    {'(a)rd', 2}, ...
    {'(r)d', 2}, ...
    {'r(d)', 2}, ...
    {'(a)b', '(a)bh', 2}, ...
    {'a(b)',  'a(b)h',  2}, ...
    {'ab(h)', 2}, ...
    {'(a)j', 3}, ...
    {'(j)jbh', '(j)jb', 3}, ...
    {'(j)bh', '(j)b', 3}, ...
    {'j(b)h', 'j(b)', 3}, ...
    {'jb(h)', 3}, ...
    };

%% ==== recode the way motif numebr and position are coded
% DONT ask me how this work.s it just does...

for i=1:length(MotifDatabase.bird)
    
    % === firsr column (which motif)
    tmp = cellfun(@(x)(x{end}), MotifDatabase.bird(i).motifs);
    
    % === second colum.. (position within mtioif)
    indtmp = find([1 diff(tmp)]);
    tmp2 = nan(size(tmp));
    tmp2(find([1 diff(tmp)])) = find([1 diff(tmp)]);
    %     tmp2(end
    
    %     counter = nan(size(tmp));
    for j=1:length(indtmp)
        
        if j==length(indtmp)
            tmp2(indtmp(j)+1:end) = tmp2(indtmp(j));
        else
            tmp2(indtmp(j)+1:indtmp(j+1)-1) = tmp2(indtmp(j));
        end
        
        
    end
    
    tmp2 = [[1:length(tmp2)] - tmp2] +1;
    out = [tmp' tmp2'];
    
    % ============= OUTPUT
    MotifDatabase.bird(i).motifs = cellfun(@(x)(x(1:end-1)), ...
        MotifDatabase.bird(i).motifs, 'UniformOutput', 0); % remove the last number
    MotifDatabase.bird(i).motifposition = out;
end


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
    motifposition_all = [];
else
    if iscell(motifregexp)
        indmotif_all = [];
        motifposition_all = [];
        for j=1:length(motifregexp)
            motifthis = motifregexp{j};
            functmp = @(X)any(strcmp(X, motifthis));
            indmotif = find(cellfun(functmp, MotifDatabase.bird(indbird).motifs));
            motifposition = MotifDatabase.bird(indbird).motifposition(indmotif, :);
            assert(length(indmotif)<2, 'PROBLEM - entered same motif in multiple entries ...');
            if isempty(indmotif)
                keyboard
                disp('PROBLEM --- this motif not entered for this bird...');
            end
            indmotif_all = [indmotif_all; indmotif];
            motifposition_all = [motifposition_all; motifposition];
        end
    else % then is a single motif
        functmp = @(X)any(strcmp(X, motifregexp));
        indmotif_all = find(cellfun(functmp, MotifDatabase.bird(indbird).motifs));
        motifposition_all = MotifDatabase.bird(indbird).motifposition(indmotif_all, :);
        assert(length(indmotif_all)<2, 'PROBLEM - entered same motif in multiple entries ...');
        if isempty(indmotif_all)
            keyboard
            disp('PROBLEM --- this motif not entered for this bird...');
        end
    end
end
