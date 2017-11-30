function DatStructByBranch = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH)
%% lt 11/29/17 - organizes data by Unique branches (defined by regexpr str)
% OUTPUTS A NEW STRUCTURE


%%
apos = 1; % assume is new analysis version, this always 1.

%% extract all branches first

ALLBRANCH = lt_neural_v2_CTXT_GetListOfBranches(ALLBRANCH);
%% PRODUCE NEW STRUCTURE ORGANIZED BY UNIQUE REGEXP STR

numbirds = length(ALLBRANCH.alignpos(apos).bird);
DatStructByBranch = struct;

for i=1:numbirds
    
    ListOfBranches = ALLBRANCH.alignpos(apos).bird(i).ListOfBranches;
    
    
    % ################################### initiate structure
    for kk=1:length(ListOfBranches)
        DatStructByBranch.bird(i).branchID(kk).regexpstr = '';
        DatStructByBranch.bird(i).branchID(kk).DAT.xdecode = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.brainregion = {};
    end
    
    
    % ################################### COLLECT DATA ND PUT INTO STRUCTURE
    numbranches = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for bb =1:numbranches
        
        numneurons = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron);
        for nn=1:numneurons
            
            % ================ collect data
            x_decode = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).xtimes;
            y_decode = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals;
            thisbranch = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).prms_regexpstr;
            
            if isempty(x_decode)
                continue
            end
            
            % ============== save
            ind = strcmp(ListOfBranches, thisbranch);
            DatStructByBranch.bird(i).branchID(ind).regexpstr = thisbranch;
            
            DatStructByBranch.bird(i).branchID(ind).DAT.xdecode = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.xdecode; ...
                x_decode];
           
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode; ...
                y_decode];
            
            location = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).NOTE_Location;
            DatStructByBranch.bird(i).branchID(ind).DAT.brainregion = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.brainregion; ...
                location];
             
        end
    end
end