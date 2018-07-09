function TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Base(TrialStruct)

%% lt 5/30/18 -convert the way baseline is indicated

%%

Numbirds = length(TrialStruct.birds);

for i=1:Numbirds
    
    Numexpt = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:Numexpt
        
        
        % ============ WHEN DID WN START?
        WNontime = TrialStruct.birds(i).exptnum(ii).BaseDays(end)+1;
        TrialStruct.birds(i).exptnum(ii).WNontime = WNontime;
        
    end
    
    % ============ REMOVE INDICATIONS OF WN DAYS
    TrialStruct.birds(i).exptnum = ...
        rmfield(TrialStruct.birds(i).exptnum, 'BaseDays');
    TrialStruct.birds(i).exptnum = ...
        rmfield(TrialStruct.birds(i).exptnum, 'WNDays');
    TrialStruct.birds(i).exptnum = ...
        rmfield(TrialStruct.birds(i).exptnum, 'WNday1');
    
end
