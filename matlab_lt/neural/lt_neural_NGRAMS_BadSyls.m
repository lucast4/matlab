function BadStrings = lt_neural_NGRAMS_BadSyls(birdname, Params)

%             if strcmp(SummaryStruct.birds(i).neurons(ii).whosbird, 'Sam');
%                BadSyls = {'0', 'X', 'x', 'i', '.'};
%             elseif strcmp(SummaryStruct.birds(i).neurons(ii).whosbird, 'Mel')
%                BadSyls = {'0', 'i', '.', 'C', 'x'};
%                if strcmp(SummaryStruct.birds(i).birdname, 'Pk35G27')
%                    BadSyls = {'0', '.', 'C', 'x'};
%                elseif strcmp(SummaryStruct.birds(i).birdname, 'O55Pu53')
%                    BadSyls = {'0', 'j', '.', 'C', 'x'};
%                end
%             end
%
%
%             wh44 -
%             '[a-z]bj', '[a-z]jn', '[a-z]bn', '[a-z]hm'


% HandCoded(1).birdname = 'B53O71';
% HandCoded(1).branches = {'[a-z]kk', '[a-z]gg', '[a-z]ca', '[a-z]aj', '[a-z]ac'};
%
% HandCoded(2).birdname = 'Pu55Pu22';
% HandCoded(2).branches = {'[a-z]ca'};
%
% HandCoded(3).birdname = 'W32Pi51';
% HandCoded(3).branches = {'[a-z]dd'};
%
% HandCoded(3).birdname = 'wh44wh39';
% HandCoded(3).branches = {'[a-z]bj', '[a-z]jn', '[a-z]bn', '[a-z]hm', 'bj[a-z]', 'cb[a-z]', 'nh[a-z]'}; % {bad labels, seq change because WN, seq change bc WN};
%
%


%% save information about ngrams to exclude because of poor labeling

if strcmp(Params.strtype, 'xaa')
    % these were checked with assumption that care about convergence points
    % ...
    % -- will remove a given ngram pair if either one contains the bad
    % string
    
    if strcmp(birdname, 'Pk35G27')
        BadStrings = {};
        
    elseif strcmp(birdname, 'B53O71')
        BadStrings = {'jkk', '[a-z]ca', '[a-z]aj', '[a-z]ac', 'i', 'gaa'};
        
    elseif strcmp(birdname, 'G26G23')
        BadStrings = {'gii', 'u'};
        
    elseif strcmp(birdname, 'G45G46')
        BadStrings = {'aae'};
        
    elseif strcmp(birdname, 'O14O15')
        BadStrings = {'aaa'};
        
    elseif strcmp(birdname, 'O55Pu53')
        BadStrings = {'j'};
        
    elseif strcmp(birdname, 'Pu55Pu22')
        BadStrings = {'cc[a-z]'};
        
    elseif strcmp(birdname, 'W15W94')
        BadStrings = {};
        
    elseif strcmp(birdname, 'W32Pi51')
        BadStrings = {};
        
    elseif strcmp(birdname, 'pu24w39')
        BadStrings = {'j', 'i'};
        
    elseif strcmp(birdname, 'pu26y2')
        BadStrings = {};
        
    elseif strcmp(birdname, 'pu44w52')
        BadStrings = {'jii'};
        
    elseif strcmp(birdname, 'W96Pi45')
        BadStrings = {'i'};
        
    elseif strcmp(birdname, 'bk7')
        BadStrings = {};
        
    elseif strcmp(birdname, 'bu77wh13')
        BadStrings = {};
        
    elseif strcmp(birdname, 'wh6pk36')
        BadStrings = {};
        
    elseif strcmp(birdname, 'br92br54')
        BadStrings = {};
        
    elseif strcmp(birdname, 'or74bk35')
        BadStrings = {};
        
    elseif strcmp(birdname, 'pu69wh78')
        BadStrings = {};
        
    elseif strcmp(birdname, 'wh44wh39')
        BadStrings = {'jjn'};
    end
end

assert(exist('BadStrings', 'var')==1, 'bird not coded ...');

% ============ add 'x', which is placeholder syl
BadStrings = [BadStrings 'x'];


