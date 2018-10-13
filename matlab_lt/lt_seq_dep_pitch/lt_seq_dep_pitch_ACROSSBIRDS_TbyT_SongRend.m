function TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_SongRend(TrialStruct)
%% lt 10/1/18 - converts datapoints to use songs instead of rends.
% will automaticalyl ceheck to make sure hasn't already been done.

if ~isfield(TrialStruct, 'ConvertedToSongAsRend')
    TrialStruct.ConvertedToSongAsRend = 0;
end

%% 

Numbirds = length(TrialStruct.birds);

%% ======= CONVERT TO USING SONG AS REND?
if TrialStruct.ConvertedToSongAsRend==0
    disp('CONVERTING FROM USING RENDS TO USING SONGS!');
    %    disp('REMOVING isWN information for now. figure out way to summarize WN hit info for song by song');
    disp('SUMMARIZING isWN information by num hits and misses');
    for i=1:Numbirds
        Numexpt = length(TrialStruct.birds(i).exptnum);
        
        for ii=1:Numexpt
            
            % ---------- SKIP IF NO DATA
            if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
                disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                continue
            end
            
            Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
            disp([num2str(i) '-' num2str(ii)]);
            % ================== go thru all syls
            for ss =1:Numsyls
                
                t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
                t_dnum = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals_datenum;
                ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
                
                if isfield(TrialStruct.birds(i).exptnum(ii).sylnum(ss), 'isCatch')
                    iscatch = TrialStruct.birds(i).exptnum(ii).sylnum(ss).isCatch;
                    if isempty(iscatch)
                        iscatch = nan(size(ff));
                    end
                else
                    iscatch = nan(size(ff));
                end
                
                if isfield(TrialStruct.birds(i).exptnum(ii).sylnum(ss), 'isWNhit')
                    isWN = TrialStruct.birds(i).exptnum(ii).sylnum(ss).isWNhit;
                else
                    isWN = nan(size(ff));
                end
                
                % ----------------- GO THRU EACH SONG AND COLLECT
                ff_songs = grpstats(ff, t, {'mean'});
                
                [~, indtmp] = unique(t);
                t_songs = t(indtmp);
                t_dnum = t_dnum(indtmp);
                
                assert(length(t_songs) == length(ff_songs));
                
                % --------- CATCH
                % assumes that catch covers entire song, so just take
                % median.
                iscatch_song = grpstats(iscatch, t, {'mean'});
                assert(all(ismember(iscatch_song(~isnan(iscatch_song)), [0 1])), 'some songs with mixed catch/nc?');
                
                % ------- WN
                % - remove trial by trial WN hit
                if (0) % DONT DO THIS, this was before I foudn way to extract hits over song.
                    if isfield(TrialStruct.birds(i).exptnum(ii).sylnum, 'isWNhit')
                        TrialStruct.birds(i).exptnum(ii).sylnum = ...
                            rmfield(TrialStruct.birds(i).exptnum(ii).sylnum, 'isWNhit');
                    end
                end
                
                % collect number of WN hits/escapes
                if all(isnan(isWN))
                    %                   wn_num = zeros(size(t_songs));
                    %                   wn_mean = zeros(size(t_songs));
                    wn_nhits = nan(size(t_songs));
                    wn_nmiss = nan(size(t_songs));
                else
                    [wn_num, wn_mean] = grpstats(isWN, t, {'numel', 'mean'});
                    wn_nhits = wn_mean.*wn_num;
                    wn_nmiss = wn_num - wn_nhits;
                    assert(all(wn_nhits == floor(wn_nhits)), 'shoud be int');
                    assert(all(wn_nmiss == floor(wn_nmiss)), 'should be int');
                end
                
                %                 if i==1 & ii==4 & ...
                %                         TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget==1
                %                     keyboard
                %                 end
                %                 if i==1 & ii==4
                %                     keyboard
                %                 end
                
                % ----------------- PUT BACK INTO STRUCTURE
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals = ff_songs;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals = t_songs;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals_datenum = t_dnum;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).isCatch = iscatch_song;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).WN_nhits = wn_nhits;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).WN_nmiss = wn_nmiss;
                
                
                % ------------- for all renditions calculate deviation from
                % recent trials
                % === method1 - fit regression line to one hour of data
                % (directly preceding this rendition...) record deviation from
                % that hour's prediction
                
            end
        end
    end
    
    TrialStruct.ConvertedToSongAsRend = 1;
end