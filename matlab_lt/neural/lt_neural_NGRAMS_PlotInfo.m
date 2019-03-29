maxbirds = max(OUTSTRUCT.All_birdnum);
for i=1:maxbirds
    
    bname = SummaryStruct.birds(i).birdname;
    
    indsthis = find(OUTSTRUCT.All_birdnum==i & ...
        OUTSTRUCT.All_diffsyl_PairType == find(strcmp(OUTSTRUCT.PairTypesInOrder, pairtype)));
    
    
    neurlist = unique(OUTSTRUCT.All_neurnum(indsthis));
    
    disp(' ========================================= ');
    disp(bname);
    
    for nn=neurlist'
        bregion = SummaryStruct.birds(i).neurons(nn).NOTE_Location;
        disp(['neur: ' num2str(nn) ' -- ' bregion]);
    end
    
    % ===== display pairs of ngrams
    ngrampairs = OUTSTRUCT.All_ngramstring_inorder(indsthis,:);
    N = OUTSTRUCT.All_N(indsthis);
    
    [inds_out, inds_unique, X_cell] = lt_tools_grp2idx({ngrampairs(:,1), ngrampairs(:,2)});
    [~, indstmp] = unique(inds_out);
    ngrampairs = ngrampairs(indstmp,:);
    N = N(indstmp);
    
    %     ngrampairs = unique(ngrampairs, 'rows');
    
    fname = [savedir '/ngramlist.txt'];
    fid = fopen(fname, 'w');
    cellfun(@(x)fwrite(fid, x), ngrampairs);
    %     fwrite(fid, ngrampairs);
    fclose(fid);
    disp(['Pairtype: ' num2str(pairtype)]);
    
    for j=1:size(ngrampairs,1)
        disp([ngrampairs(j,:) ' -- N= ' num2str(N(j))]);
    end
    
    
    % ====================== LOAD BIRD DATA
    cd(savedir)
    tmp = load(['bird' num2str(i) '.mat']);
    ngramslistall = {};
    Nall = [];
    STDmax = []; % max std of on-on [max over all syls in motif]
    for j=1:length(tmp.birdstruct.neuron)
        
        for jj=1:length(tmp.birdstruct.neuron(j).ngramlist)
            if isempty(tmp.birdstruct.neuron(j).ngramnum(jj).DAT)
                continue
            end
            %            off = [tmp.birdstruct.neuron(j).ngramnum(jj).DAT.motifsylOff];
            on = [tmp.birdstruct.neuron(j).ngramnum(jj).DAT.motifsylOn];
            assert(Params.alignsyl==2, 'i assuem that deviation from syl 2 is value to use as mean in cv..')
            tmp1 = std(on, [], 1);
            %            tmp2 = diff(mean(on,1), [], 2);
            %            tmp1([1 3:end])./tmp2
            
            % ====== save all
            ngramslistall = [ngramslistall; tmp.birdstruct.neuron(j).ngramnum(jj).regexprstr];
            Nall = [Nall; length(tmp.birdstruct.neuron(j).ngramnum(jj).DAT.tvals)];
            STDmax = [STDmax; max(tmp1)];
            
        end
    end
    [~, indsort] = sort(ngramslistall);
    ngramslistall = ngramslistall(indsort);
    Nall = Nall(indsort);
    STDmax = STDmax(indsort);
    
    Nlist = round(grpstats(Nall, ngramslistall));
    [stdmaxlist, nmot] = grpstats(STDmax, ngramslistall, {'mean', 'numel'});
    ngramlist = unique(ngramslistall);
    
    % =========== disp all motifs
    disp('All Ngrams, mean sample size, mean of max onset-onset standard deviation, Nneur with this motif: ');
    %     indsthis = find(OUTSTRUCT.All_birdnum==i);
    %     ngrams = OUTSTRUCT.All_ngramstring_inorder(indsthis,:);
    %
    %     disp(unique(ngrams(:)));
    for j=1:length(ngramlist)
        disp([ngramlist{j} ' -- N= ' num2str(Nlist(j)) ' -- on-onSTD = ' num2str(round(1000*stdmaxlist(j))) ' ms ---Nneur: ' num2str(nmot(j))]);
    end
end
