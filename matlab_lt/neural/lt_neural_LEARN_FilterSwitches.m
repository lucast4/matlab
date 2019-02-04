function indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday, birdtoget)
%% lt 12/16/18 - FIND DATA THAT MATCH CRITERIA OF SWITCH TYPE
% == will consider it a match if EVERY target motif for a given switch
% pasess criterion.

% swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
% % for a given switch did not have a previous switch on the same day
% leave expty to get all switches.

if ~exist('firstswitchofday', 'var')
    firstswitchofday = 0; % if 1, then must be first switch out of any target.
end

if ~exist('birdtoget', 'var')
    birdtoget = [];
end
%% hand entered things that should not be included...

% === wh72pk12, RALMANLearn3, switch 7 should be excluded becuase show
% weird opposite direction learning.


%%
indtoget_b_e_s = [];
for j=1:length(SwitchStruct.bird)
    birdname = SwitchStruct.bird(j).birdname;
    
    if ~isempty(birdtoget)
        if ~ismember(j, birdtoget)
            disp(['SKIPPING bird ' num2str(j)]);
            continue
        end
    end
    for jj=1:length(SwitchStruct.bird(j).exptnum)
        exptname = SwitchStruct.bird(j).exptnum(jj).exptname;
        
        for ss=1:length(SwitchStruct.bird(j).exptnum(jj).switchlist)
            
            % ######################### IS THIS AD HOC
            % EXPERIMENT TO EXCLUDE? % if so, automatically
            % exclude it.
            if strcmp(birdname, 'wh72pk12') & strcmp(exptname, 'RALMANLearn3') ...
                    & ss==7
                continue
                % excluded becuase show
                % weird opposite direction learning.
%             elseif strcmp(birdname, 'gr48bu5') & strcmp(exptname, 'RALMANLearn3')
%                 continue
%             elseif strcmp(birdname, 'wh72pk12') & strcmp(exptname, 'RALMANLearn6')
%                 continue
                % --- no learning.
            end
            
            % ##################################
            % === what other switches are on thisday?
            tmp = [SwitchStruct.bird(j).exptnum(jj).switchlist(1:ss-1).switchdnum];
            
            sw_other_today = find(floor(tmp) == ...
                floor(SwitchStruct.bird(j).exptnum(jj).switchlist(ss).switchdnum));
            
            % ===
            learningContingencies = SwitchStruct.bird(j).exptnum(jj).switchlist(ss).learningContingencies;
            
            allconting = learningContingencies(2:2:end);
            
            tmp = [];
            for k=1:length(allconting)
                
                tmp = [tmp any(cellfun(@(x)(all(x==allconting{k})), swtoget))];
                
                % --- check if this motif had a previous switch today
                if firstswitchfortarget_withinday==1
                    motifthis = learningContingencies{2*k-1};
                    for kk=sw_other_today
                        if any(strcmp(SwitchStruct.bird(j).exptnum(jj).switchlist(kk).learningContingencies(1:2:end), ...
                                motifthis))
                            % then this motif has been in a previous switch todya..
                            tmp(end) = 0;
                            break
                        end
                    end
                end
                
            end
            
            if firstswitchofday==1
                if ~isempty(sw_other_today)
                    % then this is not toda's first switch
                    continue
                end
            end
            
            if all(tmp) | isempty(swtoget)
                % then there is at least one of the swtogets for which all
                % the target syls match...
                indtoget_b_e_s = [indtoget_b_e_s; [j jj ss]];
                continue
            end
            
            
            
            
            
            %                 cellfun(@(x)ismember(swtoget, x), allconting)
            %
            %                 cellfun(@(x)ismember(x, swtoget), allconting)
            %
            
            
            %             tmp = cellfun(@(y)(all(cellfun(@(x)(all(x==y)), allconting))), swtoget);
            %
            %             if any(tmp)
            %                 % then there is at least one of the swtogets for which all
            %                 % the target syls match...
            %                 indtoget_b_e_s = [indtoget_b_e_s; [j jj ss]];
            %             end
            
            
            %             strtoplot = '';
            %             for k=1:length(learningContingencies)/2
            %                strtoplot = [strtoplot ' -- ' [learningContingencies{2*k-1} ' [' num2str(learningContingencies{2*k}) ']']];
            %             end
            %             bname = SwitchStruct.bird(j).birdname;
            %             ename = SwitchStruct.bird(j).exptnum(jj).exptname;
            %             disp(' ');
            %             disp([bname '-' ename '-sw' num2str(ss)]);
            %             disp(strtoplot);
        end
        
    end
end



