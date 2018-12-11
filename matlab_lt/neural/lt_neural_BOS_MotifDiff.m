function [AllUnitsMotifDiff, birdname] = lt_neural_BOS_MotifDiff(DATSTRUCT, SummaryBOS, PARAMS, ...
    exptlist, bostype, bregion, motif1, motif2)

%% lt 12/3/18 - for a given pair of motifs, extract across units the difference between motifs (e.g. FR)


[AllUnits1, birdname] = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS, PARAMS, ...
    exptlist, bostype, bregion, motif1);
[AllUnits2, birdname] = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS, PARAMS, ...
    exptlist, bostype, bregion, motif2);

% ========= pair up units
maxexpt = max(exptlist);
maxunit = max([AllUnits1.unitnum; AllUnits2.unitnum]);

AllUnitsMotifDiff.frdiff = {};
AllUnitsMotifDiff.t = {};

for j=exptlist
    for jj=1:maxunit
       
        % === FIND
        ind1 = AllUnits1.exptnum==j & AllUnits1.unitnum==jj;
        ind2 = AllUnits2.exptnum==j & AllUnits2.unitnum==jj;
       
        if ~any(ind1)
            continue
        end
        
        assert(sum(ind1)==1); assert(sum(ind2)==1);
        
        fr1 = mean(AllUnits1.FRmat{ind1},2);
        t1 = AllUnits1.FRmat_t{ind1};
        
        fr2 = mean(AllUnits2.FRmat{ind2},2);
        t2 = AllUnits2.FRmat_t{ind2};
        
        % === clip longer one to match shorter one
        maxinds = min([length(t1), length(t2)]);
        fr1 = fr1(1:maxinds);
        t1 = t1(1:maxinds);
        fr2 = fr2(1:maxinds);
        t2 = t2(1:maxinds);
        
        frdiff = fr1-fr2;
        assert(t1(end)==t2(end));
        
        % ============ save output
        AllUnitsMotifDiff.frdiff = [AllUnitsMotifDiff.frdiff; frdiff];
        AllUnitsMotifDiff.t = [AllUnitsMotifDiff.t; t2];
        
        %         % ==== only keep the desired unit
%         AllUnits1_SUB = lt_structure_subsample_all_fields(AllUnits1, find(ind1), 1);
%         AllUnits2_SUB = lt_structure_subsample_all_fields(AllUnits2, find(ind2), 1);
%         
%         % ==== extract data
%         stats1 = lt_neural_BOS_StatsUnits(AllUnits1_SUB, SummaryBOS, PARAMS);
%         stats2 = lt_neural_BOS_StatsUnits(AllUnits2_SUB, SummaryBOS, PARAMS);
%         
%         % === subtract
%         fr1 = stats1.frate_unitmean{1};
%         t1 = AllUnits1_SUB.FRmat_t
    end
end

