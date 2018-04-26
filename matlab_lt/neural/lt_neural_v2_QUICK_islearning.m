function [islearning, LearnSummary, switchtime] = lt_neural_v2_QUICK_islearning(birdname, exptname, extractlearninfo)
%% updated (lt) 4/23, so that have to have at least 1day of no WN 
% to be considered a WN off epoch...

%%
% extractlearninfo =1, then gets summary fo learning for this expt. if 0, then just tells you if islearning.
% [note: this only runs if actually is laernig]

if ~exist('extractlearninfo', 'var')
    extractlearninfo = 0;
end

%% lt - decides if this was a learning expt

LearnStruct = lt_neural_v2_LoadLearnMetadat;

% is this learning expt?
birdindtmp = strcmp({LearnStruct.bird.birdname}, birdname);
if ~any(birdindtmp)
    islearning = [];
    LearnSummary = [];
    switchtime = [];
    return
    
end

islearning = any(strcmp([LearnStruct.bird(birdindtmp).info(1,:)], exptname));

% --- optional: output info about this learn expt
if extractlearninfo==1
    LearnSummary = struct;
    
    if islearning==1
        
        exptinds = strcmp([LearnStruct.bird(birdindtmp).info(1,:)], exptname);
        
        LearnDat = LearnStruct.bird(birdindtmp).info(:, exptinds);
        
        % --- extract vector of datenums (switches) and statuses (pre
        % and post)
        numtargs = size(LearnDat, 2);
        
        for i =1:numtargs
            
            targsyl = LearnDat{2, i};
            LearnSummary.targnum(i).targsyl = targsyl;
            
            switches = LearnDat(3:end, i);
            
            switches = switches(~cellfun('isempty', switches));
            
            % - for each switch extract datenum and statuses
            numswitches = length(switches);
            
            for j=1:numswitches
                
                dnum = datenum(switches{j}(1:14), 'ddmmmyyyy-HHMM');
                statuspre = switches{j}(16:17);
                statuspost = switches{j}(19:20);
                
                LearnSummary.targnum(i).switches(j).datenum = dnum;
                LearnSummary.targnum(i).switches(j).statuspre = statuspre;
                LearnSummary.targnum(i).switches(j).statuspost = statuspost;
                
            end
        end
    end
end

%% ===== output time of WN onset (earliest time across all targs)
if extractlearninfo==1
    switchtime = [];
    if islearning==1
        numtargs = length(LearnSummary.targnum);
        for jj=1:numtargs
            
            
            if find(strcmp({LearnSummary.targnum(jj).switches.statuspre}, 'Of'), 1, 'first') ~=1
                % then this xperiments started with something other than WN
                % ,,, throw out all data
                indtmp = [];
            else
                
                % -- find the latest switch such that 1) the preceding
                % epoch has WN off and 2) there is at least a day of WN
                % being off.
                indtmp = strcmp({LearnSummary.targnum(jj).switches.statuspre}, 'Of') ... % switches with pre being off
                    & diff([0 LearnSummary.targnum(jj).switches.datenum])>1;  % switches preceded by >1 day of stability
                indtmp = max(find(indtmp)); % the latest switch (to take most data);
                
                % old version, runs into issue if e.g. on--> off--> on
                % during experiment, then will keep data for that off, even
                % if it is e.g. only an hour long...
                if (0)
                indtmp = find(strcmp({LearnSummary.targnum(jj).switches.statuspre}, 'Of'), 1, 'last');
                end
            end
            
            if isempty(indtmp)
                % then WN was always on
                swtimethis = 0; % make this 0 so this will always be earliest, and therefore must throw out all data
            else
                swtimethis = LearnSummary.targnum(jj).switches(indtmp).datenum;
            end
            
            
            % time when WN began
            switchtime = min([switchtime swtimethis]); % get earliest time across all targets
        end
    else
    end
end