%% ================= get bregion for each datapt
All_Bregion = cell(length(OUTSTRUCT.All_birdnum),1);
for i=1:length(OUTSTRUCT.All_birdnum)
    bnum = OUTSTRUCT.All_birdnum(i);
    neur = OUTSTRUCT.All_neurnum(i);
    disp(num2str(i));
    All_Bregion{i} = SummaryStruct.birds(bnum).neurons(neur).NOTE_Location;
end
OUTSTRUCT.All_Bregion = All_Bregion;

%% ================= get all pairwise distances during premotor window
if (0) % OLD VERSION --- this works with NGRAMSTRUCT. new version does not since
    % filesize too large.
    
    % if ever want to use this version, need to run the script in here.
    % this extracts summary arrays.
    lt_neural_NGRAMS_GetDatOldVersion;
    
end
%% ================== convert motifpair types to groups

if strcmp(Params.strtype, 'xaa')
    PairTypesInOrder = {...,
        '0  0  1', ...
        '1  0  0', ...
        '1  0  1', ...
        '0  1  0', ...
        '0  1  1', ...
        '1  1  0', ...
        '1  1  1'}; % in order to be plotted
elseif strcmp(Params.strtype, 'xaaa')
    PairTypesInOrder = {...,
        '0  0  0  1', ...
        '1  0  0  0', ...
        '1  0  0  1', ...
        '0  1  0  0', ...
        '0  1  0  1', ...
        '1  1  0  0', ...
        '1  1  0  1', ...
        '0  0  1  1', ...
        '0  0  1  0', ...
        '0  1  1  0', ...
        '0  1  1  1', ...
        '1  1  1  0', ...
        '1  0  1  0', ...
        '1  0  1  1', ...
        '1  1  1  1', ...
        }; % in order to be plotted
end
All_diffsyl_string = num2str(double(All_diffsyl_logical));
All_diffsyl_string = mat2cell(All_diffsyl_string, ones(size(All_diffsyl_logical,1),1));

[~, All_diffsyl_PairType] = ismember(All_diffsyl_string, PairTypesInOrder);
assert(all(strcmp(All_diffsyl_string, PairTypesInOrder(All_diffsyl_PairType)')), 'asdfas');

% ================= PUT INTO STRUCT
OUTSTRUCT.PairTypesInOrder = PairTypesInOrder;
OUTSTRUCT.All_diffsyl_PairType = All_diffsyl_PairType;


%% ============ [FIGURE OUT MISLABELED SYLLABLES]

ncases = length(OUTSTRUCT.All_birdnum);
OUTSTRUCT.All_BadSyls = nan(ncases,1);
for i=1:ncases
    disp(num2str(i));
    % -- birdname
    birdname = SummaryStruct.birds(OUTSTRUCT.All_birdnum(i)).birdname;
    badstrings = lt_neural_NGRAMS_BadSyls(birdname, Params);
    
    % -- if either of then grams contains any of the bad strings, then
    % mark as bad
    ngramstrings = OUTSTRUCT.All_ngramstring_inorder(i,:);
    
    %     if strcmp(birdname, 'pu26y2')
    %         keyboard
    %     end
    anybad = [];
    for j=1:2
        tmp = regexp(ngramstrings{j}, badstrings);
        anybad = any([anybad ~cellfun(@isempty, tmp)]);
    end
    
    OUTSTRUCT.All_BadSyls(i) = anybad;
end


% =========== for each bird, list all motif pairs that are bad
if (0)
    numbirds = max(OUTSTRUCT.All_birdnum);
    for j=1:numbirds
        inds = OUTSTRUCT.All_birdnum==j & OUTSTRUCT.All_BadSyls==1;
        motifpairs = OUTSTRUCT.All_ngramstring_inorder(inds,:);
        birdname = SummaryStruct.birds(j).birdname;
        
        disp([' ==================== ' birdname]);
        
        motifpairstrings = cell(size(motifpairs,1), 1);
        for jj=1:size(motifpairs,1)
            motifpairstrings{jj} = [motifpairs{jj,1} '-' motifpairs{jj,2}];
        end
        disp(unique(motifpairstrings));
        disp(' ---- ignoring any with "x"');
        tmp = unique(motifpairstrings);
        indsnox = cellfun(@isempty, regexp(tmp, 'x'));
        disp(tmp(indsnox));
    end
    
end

%% ============ REMOVE EMPTY FIELDS FROM OUTSTURCT

OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT);
%% ============ REMOVE BAD LABEL PAIRS FROM DATASET

% ==== save original outsturct
OUTSTRUCT_orig = OUTSTRUCT;

% ==== save this for later
PairTypesInOrder = OUTSTRUCT.PairTypesInOrder;

% ==== remove this field, since is smaller vector.
OUTSTRUCT = rmfield(OUTSTRUCT, 'PairTypesInOrder');

% ==== go thru all fields and only keep the good syls
indstokeep = ~OUTSTRUCT.All_BadSyls;
fnames = fieldnames(OUTSTRUCT);

for j=1:length(fnames)
    
    ytmp = OUTSTRUCT.(fnames{j});
    
    try
        ytmp = ytmp(indstokeep, :);
        OUTSTRUCT.(fnames{j}) = ytmp;
    catch err
        disp('ok');
        try ytmp = ytmp(:, indstokeep);
            OUTSTRUCT.(fnames{j}) = ytmp';
        catch err
            disp('why error?');
        end
    end
end

% === put pairtypes back in
OUTSTRUCT.PairTypesInOrder = PairTypesInOrder;