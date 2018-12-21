function DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_UNDIR, ListOfDirs_DIR, ...
    ListOfBatch, MotifsToExtract, gethitsyls)

if ~exist('gethitsyls', 'var')
    gethitsyls=1;
end

%% lt 11/14/17 - extracts results from lt_batchsong_calcFF


%% ============ cutoffs for determining what is a hit.

durmin = 0.02; % duration of wn, sec
durmax = 0.1;
npeaksmin = 100;
npeaksmax = 2000;


%%
ListOfDirs_ALL = [ListOfDirs_UNDIR ListOfDirs_DIR];

%%
if (0) % old version, too slow since loads each song file once for each syl.
    DATSTRUCT = struct;
    for mm = 1:length(MotifsToExtract)
        
        regexpr_str = MotifsToExtract{mm};
        DATSTRUCT.motif(mm).motif = regexpr_str;
        
        count=1;
        for i=1:length(ListOfBatch)
            
            dirname = ListOfDirs_ALL{i};
            batchf = ListOfBatch{i};
            
            % ==== is this DIR or UNDIR song?
            if (1)
                % new version, first dirs and undir, second directoreis and DIR
                % song.
                if i <= length(ListOfDirs_UNDIR)
                    isDIR = 0;
                else
                    isDIR=1;
                end
            else
                if any(ismember(ListOfDirs_UNDIR, dirname))
                    isDIR = 0;
                elseif any(ismember(ListOfDirs_DIR, dirname))
                    isDIR = 1;
                else
                    disp('PROBLEM!!! - dir or undir?');
                end
            end
            
            cd(dirname);
            
            % --- go thru all songs in batchf
            fid = fopen(batchf);
            fname = fgetl(fid);
            
            while ischar(fname)
                disp(fname);
                
                % ================ check if calFF.mat exists
                if ~exist([fname '.calcff.mat'], 'file')
                    fname = fgetl(fid);
                    disp('==== SKIP - no notmat')
                    continue
                end
                
                % ============= extract FF and time for all syls
                notmat = load([fname '.not.mat']);
                calcff = load([fname '.calcff.mat']);
                assert(length(notmat.labels) == length(calcff.FFstruct.FFall));
                
                
                
                %% ============== GET DATA ON HITS/ESCAPES
                [~, fn, fe] = fileparts(fname);
                if strcmp(fe, '.cbin')
                    % then feedback info is in rec file
                    
                    rd=readrecf(fname);
                    
                    trigtimes = rd.ttimes; % in msec
                    iscatch = rd.iscatch;
                    
                    
                elseif strcmp(fe, '.rhd')
                    
                    % then feecbak must be extracted
                    wntime = load([fname '.wntime.mat']);
                    
                    trigtimes = wntime.wnstruct.WNonsets*1000; % convert to ms.
                    iscatch = 0; % for rhd assume there is no catch. this is true for all
                    % experiments up to now (8/2018).
                    
                    % ===================== filter to determine waht is a true
                    % hit
                    durtmp = ([wntime.wnstruct.WNoffsets] - [wntime.wnstruct.WNonsets]);
                    indstmp = durtmp>durmin & durtmp<durmax & ...
                        wntime.wnstruct.Nmaxes_withindur>npeaksmin & ...
                        wntime.wnstruct.Nmaxes_withindur<npeaksmax;
                    trigtimes = trigtimes(indstmp);
                end
                
                
                % ================== GET POSITIONS OF HIT SYLS
                hitInds = [];
                for j=1:length(trigtimes)
                    hs = find(trigtimes(j)>notmat.onsets ...
                        & trigtimes(j)<notmat.offsets);
                    
                    hitInds = [hitInds; hs];
                    %                     notmat.labels(hitsyl-2:hitsyl+2)
                end
                
                %%
                % ---------------- extract time of song
                [~, ~, fext] = fileparts(fname);
                if strcmp(fext, '.rhd')
                    [dtnum datestring]=lt_neural_fn2datenum(fname);
                elseif strcmp(fext, '.cbin')
                    dtnum = fn2datenum_eftafv4_lt(fname);
                    datestring = datestr(dtnum, 'yymmdd_HHMMSS');
                end
                
                
                % =============== find all instances of this motif in this
                % songfile
                [tokenExtents, startinds, endinds, matchlabs] = lt_batchsong_regexp(notmat.labels, regexpr_str);           %%
                if isempty(tokenExtents)
                    fname = fgetl(fid);
                    continue
                end
                assert(size(tokenExtents, 1)==1, 'needs to be horizontal');
                
                for j=tokenExtents
                    
                    % ============= extract data for this rend
                    DATSTRUCT.motif(mm).rendnum(count).ff = calcff.FFstruct.FFall(j);
                    DATSTRUCT.motif(mm).rendnum(count).syl = notmat.labels(j);
                    DATSTRUCT.motif(mm).rendnum(count).fname = fname;
                    DATSTRUCT.motif(mm).rendnum(count).dirname = dirname;
                    DATSTRUCT.motif(mm).rendnum(count).datenum_song_SecRes = dtnum;
                    DATSTRUCT.motif(mm).rendnum(count).datestr = datestring;
                    DATSTRUCT.motif(mm).rendnum(count).time_withinsong = notmat.onsets(j)/1000;
                    DATSTRUCT.motif(mm).rendnum(count).isDIR = isDIR;
                    
                    % ----------- is this a hit?
                    DATSTRUCT.motif(mm).rendnum(count).isWNhit = ismember(j, hitInds);
                    DATSTRUCT.motif(mm).rendnum(count).isCatchsong = iscatch;
                    count = count+1;
                end
                
                
                fname = fgetl(fid);
            end
        end
        
    end
    
else
    
    %%
    DATSTRUCT = struct;
    for i=1:length(ListOfBatch)
        
        % ==== is this DIR or UNDIR song?
        if (1)
            % new version, first dirs and undir, second directoreis and DIR
            % song.
            if i <= length(ListOfDirs_UNDIR)
                isDIR = 0;
            else
                isDIR=1;
            end
        else
            if any(ismember(ListOfDirs_UNDIR, dirname))
                isDIR = 0;
            elseif any(ismember(ListOfDirs_DIR, dirname))
                isDIR = 1;
            else
                disp('PROBLEM!!! - dir or undir?');
            end
        end
        
        
        % --- go thru all songs in batchf
        dirname = ListOfDirs_ALL{i};
        batchf = ListOfBatch{i};
        cd(dirname);
        fid = fopen(batchf);
        fname = fgetl(fid);
        
        while ischar(fname)
            disp(fname);
            
            % ================ check if calFF.mat exists
            if ~exist([fname '.calcff.mat'], 'file')
                fname = fgetl(fid);
                disp('==== SKIP - no notmat')
                continue
            end
            
            
            
            % ============= extract FF and time for all syls
            notmat = load([fname '.not.mat']);
            calcff = load([fname '.calcff.mat']);
            assert(length(notmat.labels) == length(calcff.FFstruct.FFall));
            
            
            
            %% ============== GET DATA ON HITS/ESCAPES
            if gethitsyls==1
                [~, ~, fe] = fileparts(fname);
                if strcmp(fe, '.cbin')
                    % then feedback info is in rec file
                    
                    rd=readrecf(fname);
                    trigtimes = rd.ttimes; % in msec
                    iscatch = rd.iscatch;
                    
                elseif strcmp(fe, '.rhd')
                    
                    % then feecbak must be extracted
                    wntime = load([fname '.wntime.mat']);
                    
                    trigtimes = wntime.wnstruct.WNonsets*1000; % convert to ms.
                    iscatch = 0; % for rhd assume there is no catch. this is true for all
                    % experiments up to now (8/2018).
                    
                    % ===================== filter to determine waht is a true
                    % hit
                    durtmp = ([wntime.wnstruct.WNoffsets] - [wntime.wnstruct.WNonsets]);
                    indstmp = durtmp>durmin & durtmp<durmax & ...
                        wntime.wnstruct.Nmaxes_withindur>npeaksmin & ...
                        wntime.wnstruct.Nmaxes_withindur<npeaksmax;
                    trigtimes = trigtimes(indstmp);
                end
                
                
                
                % ================== GET POSITIONS OF HIT SYLS
                hitInds = [];
                for j=1:length(trigtimes)
                    hs = find(trigtimes(j)>notmat.onsets ...
                        & trigtimes(j)<notmat.offsets);
                    hitInds = [hitInds; hs];
                end
            end
            
            %%
            % ===================== extract time of song
            [~, ~, fext] = fileparts(fname);
            if strcmp(fext, '.rhd')
                [dtnum, datestring]=lt_neural_fn2datenum(fname);
            elseif strcmp(fext, '.cbin')
                dtnum = fn2datenum_eftafv4_lt(fname);
                datestring = datestr(dtnum, 'yymmdd_HHMMSS');
            end
            
            
            % ========================= GO THRU ALL MOTIFS.
            for mm = 1:length(MotifsToExtract)
                
                regexpr_str = MotifsToExtract{mm};
                if i==1
                    DATSTRUCT.motif(mm).motif = regexpr_str;
                end
                
                
                % =============== find all instances of this motif in this
                % songfile
                [tokenExtents, startinds, endinds, matchlabs] = lt_batchsong_regexp(notmat.labels, regexpr_str);           %%
                if isempty(tokenExtents)
                    continue
                end
                assert(size(tokenExtents, 1)==1, 'needs to be horizontal');
                
                for j=tokenExtents
                    
                    if ~isfield(DATSTRUCT.motif(mm), 'rendnum')
                        count =1;
                    else
                        count = length(DATSTRUCT.motif(mm).rendnum)+1;
                    end
                    
                    % ============= extract data for this rend
                    DATSTRUCT.motif(mm).rendnum(count).ff = calcff.FFstruct.FFall(j);
                    DATSTRUCT.motif(mm).rendnum(count).syl = notmat.labels(j);
                    DATSTRUCT.motif(mm).rendnum(count).fname = fname;
                    DATSTRUCT.motif(mm).rendnum(count).dirname = dirname;
                    DATSTRUCT.motif(mm).rendnum(count).datenum_song_SecRes = dtnum;
                    DATSTRUCT.motif(mm).rendnum(count).datestr = datestring;
                    DATSTRUCT.motif(mm).rendnum(count).time_withinsong = notmat.onsets(j)/1000;
                    DATSTRUCT.motif(mm).rendnum(count).isDIR = isDIR;
                    
                    % ----------- is this a hit?
                    if gethitsyls==1
                        DATSTRUCT.motif(mm).rendnum(count).isWNhit = ismember(j, hitInds);
                        DATSTRUCT.motif(mm).rendnum(count).isCatchsong = iscatch;
                    end
                end
            end
            fname = fgetl(fid);
        end
        
    end
end