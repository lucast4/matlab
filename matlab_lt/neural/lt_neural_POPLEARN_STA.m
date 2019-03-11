function sta_all = lt_neural_POPLEARN_STA(spkdat, lfpdat, trialstoget, sta_wind, ...
    t_lfp)
%% lt 2/25/19 - gets all cases of spikes, spike-triggered LFP.

%     spkdat = DATSTRUCT_SPK.spike_dat{i};
%     lfpdat = DATSTRUCT_LFP.LFP_dat{DATSTRUCT_SPK.LFP_indthis(i)};
%     assert(length(spkdat)==size(lfpdat,2));
%
%     trialstoget = DATSTRUCT_SPK.inds_base_epoch{i};
% sta_wind = [-0.05 0.05]; % relative to spike, in sec % will only get spikes that are within data...

%% get spike triggered average of lfp
sta_all = [];
for j=trialstoget
    spkthis =spkdat{j};
    lfpthis = lfpdat(:,j);
    
    % ---- for each spike, get adjacent lfp data.
    nspks = length(spkthis);
    lengthtoget = (sta_wind(2)-sta_wind(1))*(1500)-1;
    for nn=1:nspks
        spktime = spkthis(nn);
        windthis = spktime + sta_wind;
        
        if windthis(1)<t_lfp(1) | windthis(2)>t_lfp(end)
            continue
        end
        
        xtoget = find(t_lfp>=windthis(1) & t_lfp<=windthis(2));
        xtoget = xtoget(1:lengthtoget);
        
        %            assert((xtoget)==lengthtoget);
        sta_all = [sta_all; lfpthis(xtoget)'];
    end
    %         sta_trials{j} = sta_tmp;
end
