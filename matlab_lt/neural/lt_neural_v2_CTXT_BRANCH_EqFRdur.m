function ALLBRANCH = lt_neural_v2_CTXT_BRANCH_EqFRdur(ALLBRANCH)

%% lt 4/24/18 - makes sure FR at all classes have same duration (rounding errors...)

%%
apos = 1;

%%
numbirds = length(ALLBRANCH.alignpos(apos).bird);
for i=1:numbirds
   numbranch = length(ALLBRANCH.alignpos(apos).bird(i).branch);
   for ii=1:numbranch
      numneuron = length(ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron);
      for nn=1:numneuron
          if isempty(ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn).FR)
              continue
          end
          disp([num2str([i ii nn])]);   
          
          numclass = length(ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn).FR.classnum);
          
          DAT = ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn);
          alldurs = [];
          % ----- get min dur
          for cc = 1:numclass
            dur = size(DAT.FR.classnum(cc).FRsmooth_rate_CommonTrialDur,1);
            if dur==0
                continue
            end
            alldurs = [alldurs dur];
          end
          
          mindur = min(alldurs);
          
          
          % =============== go in and whittle down all durs
          frfield = 'FR';
          DAT = ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn).(frfield);
          numclass = length(DAT.classnum);
          for cc =1:numclass
              if isempty(DAT.classnum(cc).FRsmooth_rate_CommonTrialDur)
                  continue
              end
             DAT.classnum(cc).FRsmooth_rate_CommonTrialDur = ...
                 DAT.classnum(cc).FRsmooth_rate_CommonTrialDur(1:mindur,:);
             
             DAT.classnum(cc).FRsmooth_xbin_CommonTrialDur = ...
                 DAT.classnum(cc).FRsmooth_xbin_CommonTrialDur(1:mindur);
          end         
          % --- put back..
          ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn).(frfield) =DAT;
          
          
          frfield = 'FR_POSCONTR';
          DAT = ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn).(frfield);
          numclass = length(DAT.classnum);
          for cc =1:numclass
              if isempty(DAT.classnum(cc).FRsmooth_rate_CommonTrialDur)
                  continue
              end
             DAT.classnum(cc).FRsmooth_rate_CommonTrialDur = ...
                 DAT.classnum(cc).FRsmooth_rate_CommonTrialDur(1:mindur,:);
             
             DAT.classnum(cc).FRsmooth_xbin_CommonTrialDur = ...
                 DAT.classnum(cc).FRsmooth_xbin_CommonTrialDur(1:mindur);
          end         
          % --- put back..
          ALLBRANCH.alignpos(apos).bird(i).branch(ii).neuron(nn).(frfield) =DAT;
          
          
      end
   end
end