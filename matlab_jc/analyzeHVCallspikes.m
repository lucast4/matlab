function [outvars]=analyzeHVCallspikes(neuron,clusters,wins)





for mls=1:length(wins);
    winstart=wins(mls);
    winstop=wins(mls)+50;

%%%%% SYLLABLE A
                % Analyze - this loop returns:
                    % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
                    % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
                % for the first syllable, the harmonic part goes from 45-70ms (includes 16ms adjustment - so really 29-54ms)
                pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
                harmstart=45;
                harmend=70;
                neuroncount=0;
                clear spikecount pitchmean 
                for thisfolder=1:11 % for the neuron/recording wire
                    goodclust=clusters{thisfolder*2};
                    for thissyllable=1 % for the syllable
                        for thisclust=1:length(goodclust)     % for the cluster
                            neuroncount=neuroncount+1; % each unit
                            count=0; % each rendition with this unit
                            for thissong=1:length(neuron(thisfolder).song)
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                    timewindow=harmstart+winstart:harmend+winstop;
                                    spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                        -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                    spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end

                % Only analyze results for the neurons that actually spike in the vicinity
                % of the harmonic stack note - the variable 'doesspike' evaluates this
                clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
                    end

                    % Calculate correlation coefficients and corresponding t-stats
                    clear rval tval
                        for neuralunit=1:length(spikecount) % for each 'neuron'
                            if doesspike(neuralunit)
                                m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                                r=m(2);
                                rval(neuralunit)=r;
                                tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                            end
                        end
                        rval(doesspike);
                        tval(doesspike);
                SyllA.rval=rval;
                SyllA.tval=tval;
                SyllA.spikecount=spikecount;
                SyllA.pitchmean=pitchmean;
                
%%%%% SYLLABLE B
                % Analyze - this loop returns:
                    % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
                    % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
                % for the second syllable, the harmonic part goes from 110-140ms (includes 16ms adjustment - so really 94-124ms)
                clear spikecount pitchmean
                pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
                harmstart=80; % quite a bit of leeway here
                harmend=160; % quite a bit of leeway here
                neuroncount=0;
                for thisfolder=1:11 % for the neuron/recording wire
                    goodclust=clusters{thisfolder*2};
                    for thissyllable=1 % for the syllable
                        for thisclust=1:length(goodclust)     % for the cluster
                            neuroncount=neuroncount+1; % each unit
                            count=0; % each rendition with this unit
                            for thissong=1:length(neuron(thisfolder).song)
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                    timewindow=harmstart+winstart:harmend+winstop;
                                    spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                        -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                    spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end

                % Only analyze results for the neurons that actually spike in the vicinity
                % of the harmonic stack note - the variable 'doesspike' evaluates this
                clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
                    end
                    clear rval tval
                    % Calculate correlation coefficients and corresponding t-stats
                        for neuralunit=1:length(spikecount) % for each 'neuron'
                            if doesspike(neuralunit)
                                m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                                r=m(2);
                                rval(neuralunit)=r;
                                tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                            end
                        end
                        rval(doesspike);
                        tval(doesspike);
                SyllB.rval=rval;
                SyllB.tval=tval;
                SyllB.spikecount=spikecount;
                SyllB.pitchmean=pitchmean;
    
    
 %%%%% SYLLABLE C
                % Analyze - this loop returns:
                    % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
                    % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
                % for the first syllable, the harmonic part goes from 45-70ms (includes 16ms adjustment - so really 29-54ms)
                pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
                harmstart=90;
                harmend=130;
                neuroncount=0;
                clear spikecount pitchmean     
                for thisfolder=1:11 % for the neuron/recording wire
                    goodclust=clusters{thisfolder*2};
                    for thissyllable=1 % for the syllable
                        for thisclust=1:length(goodclust)     % for the cluster
                            neuroncount=neuroncount+1; % each unit
                            count=0; % each rendition with this unit
                            for thissong=1:length(neuron(thisfolder).song)
                                for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                                    timewindow=harmstart+winstart:harmend+winstop;
                                    spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                        -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                                    spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end

                % Only analyze results for the neurons that actually spike in the vicinity
                % of the harmonic stack note - the variable 'doesspike' evaluates this
                clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
                    end

                    % Calculate correlation coefficients and corresponding t-stats
                    clear rval tval
                        for neuralunit=1:length(spikecount) % for each 'neuron'
                            if doesspike(neuralunit)
                                m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                                r=m(2);
                                rval(neuralunit)=r;
                                tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                            end
                        end
                        rval(doesspike);
                        tval(doesspike);
                SyllC.rval=rval;
                SyllC.tval=tval;
                SyllC.spikecount=spikecount;
                SyllC.pitchmean=pitchmean;
   
mnabsA(mls)=median(abs(SyllA.tval(SyllA.tval~=0 & ~isnan(SyllA.tval))));
mnabsB(mls)=median(abs(SyllB.tval(SyllB.tval~=0 & ~isnan(SyllB.tval))));
mnabsC(mls)=median(abs(SyllC.tval(SyllC.tval~=0 & ~isnan(SyllC.tval))));
outvars.mnA(mls)=median((SyllA.tval(SyllA.tval~=0 & ~isnan(SyllA.tval))));
outvars.mnB(mls)=median((SyllB.tval(SyllB.tval~=0 & ~isnan(SyllB.tval))));
outvars.mnC(mls)=median((SyllC.tval(SyllC.tval~=0 & ~isnan(SyllC.tval))));
outvars.seA(mls)=std((SyllA.tval(SyllA.tval~=0 & ~isnan(SyllA.tval))))/sqrt(sum(SyllA.tval~=0 & ~isnan(SyllA.tval)));
outvars.seB(mls)=std((SyllB.tval(SyllB.tval~=0 & ~isnan(SyllB.tval))))/sqrt(sum(SyllB.tval~=0 & ~isnan(SyllB.tval)));
outvars.seC(mls)=std((SyllC.tval(SyllC.tval~=0 & ~isnan(SyllC.tval))))/sqrt(sum(SyllC.tval~=0 & ~isnan(SyllC.tval)));
outvars.neuroncountA(mls)=sum(SyllA.tval~=0 & ~isnan(SyllA.tval));
outvars.neuroncountB(mls)=sum(SyllB.tval~=0 & ~isnan(SyllB.tval));
outvars.neuroncountC(mls)=sum(SyllC.tval~=0 & ~isnan(SyllC.tval));
end
