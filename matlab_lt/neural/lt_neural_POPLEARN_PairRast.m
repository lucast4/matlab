BirdExptPairsToPlot = {};
SwitchToPlot = [2];
TypeOfPairToPlot = {'LMAN-RA'}; % e.g. 'LMAN-RA' (in alphabetical order)

numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        % ----------------- ONLY PLOT SPECIFIC BIRD?
        if ~isempty(BirdExptPairsToPlot)
           
            ind1 = find(strcmp(BirdExptPairsToPlot, birdname));
            ind2 = find(strcmp(BirdExptPairsToPlot, exptname));
            
            if ~any(ind1+1 == ind2)
                disp(['SKIPPED ' birdname '-' exptname]);
                continue
            end
            
        end
        
        % ----------------- GO THRU ALL SWITCHES
        for iii=1:numswitches
            
            if ~isempty(SwitchToPlot)
               if ~any(SwitchToPlot == iii)
                   continue
               end
            end
            
            % ---- for this switch, figure out which populations have data
            % overlapping the onset (i.e. has data both pre and post swictch)
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
            numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
            
            for ss = 1:numsets
                songfiles = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_songfiles{ss};
                songtimes = datenum(songfiles, 'yymmdd_HHMMSS');
                
                inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                
                if isempty(inds_pre) | isempty(inds_post)
                    continue
                else
                    disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
                end
                
                
                
                % ############################################### ANALYSIS/PLOTS
                DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss);
                motiflist = {DAT.motif.regexpstr};
                neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
                
                % ============ for each pair of neurons, plot paired
                % rasters
                % -- go thru all pairs of neurons, only plot if is desired
                % type of pair
                for j=1:length(neurlist)
                    for jj=j+1:length(neurlist)
                   
                        n1 = neurlist(j);
                        n2 = neurlist(jj);
                        
                        % ----- check what pair of brain region
                        
                        
                    end
                end
                
                
            end
        end
    end
end
