function tokList = lt_neural_QUICK_GetTokens(strIn, posIn)
%% lt 5/25/19 - given a string and index locations, extracts all potential tokens (up to length 2)

% e.g input:
% strIn = 'abc'; posIn=1:2

% output:
% tokList =
% 
%   7×1 cell array
% 
%     'a'
%     '(a)'
%     '(a)b'
%     'b'
%     '(b)'
%     'a(b)'
%     '(b)c'

%%
assert(size(posIn,1)==1) % must be 1 x N

tokList = {};
for i=posIn
    
    % ===== 1) length 1
    tokList = [tokList; strIn(i)];
%     tokList = [tokList; '(' strIn(i) ')'];
    
    % ===== 2) length 2 (preceding and following)
    if i>1
    tokList = [tokList; [strIn(i-1) '(' strIn(i) ')']];
    end
    if i<length(strIn)
        tokList = [tokList; ['(' strIn(i) ')' strIn(i+1)]];
    end
            
    
end