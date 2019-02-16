%% LT 1/13/15 - v3: simply writes batch files without queruying user for anything.
% LT modified 1/26/15 to also get out songs with catch FB. all included in
% rec file.
% Also modified to categorize FB songs into catchs songs and notcatch

function lt_rec_files_find_FB_v3(batch)
%% INPUTS:
% batch='batch.catch'; % batch name to process

%%

directory=pwd;

% print rec files into one batch file
rec_files.batch=batch;

% Initiate files
frec_FB=fopen([batch '.rec_FB'],'w'); % only those with feedback give (i.e. template match)
frec_noFB=fopen([batch '.rec_noFB'],'w'); % only those with no feedback
frec_FB_MoreThanOne=fopen([batch '.rec_FB_MoreThanOne'],'w');
% frec_FB_OnlyOne=fopen([batch '.rec_FB_OnlyOne'],'w');
frec_FB_catch=fopen([batch '.rec_FB.catch'],'w'); % only those with feedback give (i.e. template match)
frec_FB_notcatch=fopen([batch '.rec_FB.notcatch'],'w'); % only those with feedback give (i.e. template match)




rec_files.DateNum_songs_with_FB=[];
rec_files.DateNum_songs_without_FB=[];
% rec_files.DateNum_songs_with_FB_MoreThanOne=[];
% rec_files.DateNum_songs_with_FB_OnlyOne=[];




%% get rec files from batch, including time info.

% fid=fopen(batch,'r');
% clock=1;
% 
% while 1
%     fn=fgetl(fid);
%     if ~ischar(fn),
%         break;
%     end
%     suffix_start_index=findstr(fn,'.cbin');
%     rec_name_from_cbin=[fn(1:suffix_start_index-1) '.rec'];
%     rec_files_list(clock)=dir(rec_name_from_cbin);
%     clock=clock+1;
% end



fid=fopen(batch,'r');
fn = fgetl(fid);
rec_files_list = [];
while ischar(fn)
    suffix_start_index=findstr(fn,'.cbin');
    rec_name_from_cbin=[fn(1:suffix_start_index-1) '.rec'];
    
    rec_files_list = [rec_files_list dir(rec_name_from_cbin)];
%     
%     rec_files_list(clock)=dir(rec_name_from_cbin);
    fn = fgetl(fid);
end


%%


for i=1:length(rec_files_list);
    frecfile{i}=fileread(rec_files_list(i).name);
    tmp1=findstr(frecfile{i},'FB'); % actual FB
    tmp2=findstr(frecfile{i},'catch #'); % catch FB
    FB_indices{i}=sort([tmp1 tmp2]);
    
    % is this catch song?
    tmp3=findstr(frecfile{i},'Catch Song'); 
    CatchSong=str2num(frecfile{i}(tmp3+13)); % 1 is catch song, 0 is not
    
    
    % convert name of rec file to end in .cbin (so can use batch with    % evsonganaly
    fn=rec_files_list(i).name;
    suffix_start_index=findstr(fn,'.rec');
    rec_name_converted_cbin=[fn(1:suffix_start_index-1) '.cbin'];
    
    if length(FB_indices{i})>0;
        fprintf(frec_FB,'%s\n',rec_name_converted_cbin);
        rec_files.DateNum_songs_with_FB=[rec_files.DateNum_songs_with_FB; rec_files_list(i).datenum];
        
        if CatchSong==1; % then is catch song;
                   fprintf(frec_FB_catch,'%s\n',rec_name_converted_cbin);
        elseif CatchSong==0;
            fprintf(frec_FB_notcatch,'%s\n',rec_name_converted_cbin);
        end
    end
    
%     if length(FB_indices{i})>1
%         fprintf(frec_multFB,'%s\n',rec_name_converted_cbin);
%     end
    if length (FB_indices{i})>1;
        fprintf(frec_FB_MoreThanOne,'%s\n',rec_name_converted_cbin);
%         rec_files.DateNum_songs_with_FB_MoreThanOne=[rec_files.DateNum_songs_with_FB_MoreThanOne; rec_files_list(i).datenum];
    end
    if length(FB_indices{i})==0;
        fprintf(frec_noFB,'%s\n',rec_name_converted_cbin);
        rec_files.DateNum_songs_without_FB=[rec_files.DateNum_songs_without_FB; rec_files_list(i).datenum];
        
    end
%     if length(FB_indices{i})==1;
%         fprintf(frec_FB_OnlyOne,'%s\n',rec_name_converted_cbin);
%         rec_files.DateNum_songs_with_FB_OnlyOne=[rec_files.DateNum_songs_with_FB_OnlyOne; rec_files_list(i).datenum];
%     end
    
    




end

rec_files.FB_indices=FB_indices;

fclose(frec_FB);
fclose(frec_noFB);
fclose(frec_FB_MoreThanOne);
% fclose(frec_FB_OnlyOne);

disp(['Num songs with FB: ' num2str(length(rec_files.DateNum_songs_with_FB))]);

%% Getting times of all FB rec files
clear feedback_FB_times_HHMMSS

for i=1:length(rec_files.DateNum_songs_with_FB);
    rec_files.feedback_FB_times_HHMMSS=datestr(rec_files.DateNum_songs_with_FB,'HHMMSS');
    rec_files.feedback_FB_times_hours=str2num(rec_files.feedback_FB_times_HHMMSS(:,1:2))+str2num(rec_files.feedback_FB_times_HHMMSS(:,3:4)).*(1/60) + str2num(rec_files.feedback_FB_times_HHMMSS(:,5:6)).*(1/3600);
end


%% plot singing times (i.e. FB times)
close all
figure(1); hold on;
plot(rec_files.feedback_FB_times_hours, 1,'og');
% for i=1:num_of_episodes;
%     line(rec_files.feeding_times_scored_hours_SyncedToEvtafTime(i,:), [0.98 0.98]);
% end

title('times when .rec was made and included evtaf FB (i.e. either real song (if FB is contingent on song) or manual trigger');
ylim([0.95 1.03]);

% put number of songs into a structure
rec_files.amount_songs_any_FB.total_songs=length(rec_files.feedback_FB_times_hours);
rec_files.amount_songs_any_FB.songs_per_hour=length(rec_files.feedback_FB_times_hours)/(max(rec_files.feedback_FB_times_hours)-min(rec_files.feedback_FB_times_hours));

%% take only the real songs (i.e. excluding FB songs that are actually just manual FB triggers)
% OIBSOLETE

if (0)
    question=input('sorting out real songs from manual triggered FB files: do you think all real songs had more than one hit of FB? (while most manual triggers had only one? if so, enter y (no apostrophe). Otherwise: n ','s');
    if question=='y';
        disp('FB instances > 1 implies .rec is a real song');
        
        clear feedback_song_times_HHMMSS
        
        for i=1:length(rec_files.DateNum_songs_with_FB_MoreThanOne);
            rec_files.feedback_song_times_HHMMSS=datestr(rec_files.DateNum_songs_with_FB_MoreThanOne,'HHMMSS');
            rec_files.feedback_song_times_hours=str2num(rec_files.feedback_song_times_HHMMSS(:,1:2))+str2num(rec_files.feedback_song_times_HHMMSS(:,3:4)).*(1/60) + str2num(rec_files.feedback_song_times_HHMMSS(:,5:6)).*(1/3600);
        end
        
        %plot
        figure(2); hold on;
        plot(rec_files.feedback_song_times_hours, 1,'og');
        
        title('times when song occured (i.e. .rec file with >1 FB (so less likely to be manual trigger) was recorded');
        
        
        % put number of songs into a structure
        rec_files.amount_songs_MoreThanOneFB.total_songs=length(rec_files.feedback_song_times_hours)
        rec_files.amount_songs_MoreThanOneFB.songs_per_hour=length(rec_files.feedback_song_times_hours)/(max(rec_files.feedback_song_times_hours)-min(rec_files.feedback_song_times_hours))
        
    elseif question=='n'
        
        disp('find a way to plot only songs')
        
    end
    
end



%% EATING - obsolete.

if (0)
    if input('do you want to also plot eating episodes (y or n)? ','s') == 'y';
        values=[];
        
        yes_or_no=input('if you have not already, make a matrix storing the start and stop times for the bird pecking at the food.  This is manually scored using video.  matrix is n x 2, where n= number of feeding episodes.  entries are 6 digits: HHMMSS, in video time.  later you will enter conversion (video to evtaf time).  save this matrix as "feeding_times_scored.mat" in the day folder.  type y if you want to start entering times, or n if you have already save the matrix','s');
        rec_files.lag_time=input('what is the lag between video and evtaf? (i.e. how much time do you subtract from the video time (start=0:00) to get the evtaf time? (enter as HHMMSS)','s');
        
        
        if yes_or_no=='y';
            num_of_episodes=input('how many episodes? ');
            rec_files.feeding_times_scored=[];
            for i=1:num_of_episodes;
                values(i,:)=input(['enter episode ' num2str(i) ' values (e.g. [102657 102745])']);
                rec_files.feeding_times_scored(i,:)=values(i,:);
            end
            rec_files.feeding_times_scored=num2str(rec_files.feeding_times_scored);
        end
        if yes_or_no=='n';
            load('feeding_times_scored.mat')
            rec_files.feeding_times_scored=num2str(feeding_times_scored);
        end
        rec_files.feeding_times_scored_hours(:,1)=str2num(rec_files.feeding_times_scored(:,1:2))+str2num(rec_files.feeding_times_scored(:,3:4)).*(1/60) + str2num(rec_files.feeding_times_scored(:,5:6)).*(1/3600)
        rec_files.feeding_times_scored_hours(:,2)=str2num(rec_files.feeding_times_scored(:,9:10))+str2num(rec_files.feeding_times_scored(:,11:12)).*(1/60) + str2num(rec_files.feeding_times_scored(:,13:14)).*(1/3600)
        
        
        % subtract lag times from the feeding times get food and song in sync.
        rec_files.lag_time_hours=str2num(rec_files.lag_time(:,1:2))+str2num(rec_files.lag_time(:,3:4)).*(1/60) + str2num(rec_files.lag_time(:,5:6)).*(1/3600)
        rec_files.feeding_times_scored_hours_SyncedToEvtafTime=rec_files.feeding_times_scored_hours-rec_files.lag_time_hours;
        
        % Add feeding times to the plot of FB times (i.e. FB times)
        figure(1);
        for i=1:num_of_episodes;
            line(rec_files.feeding_times_scored_hours_SyncedToEvtafTime(i,:), [0.98 0.98]);
        end
        
        
        
        % Add feeding times to the plot of singing times (i.e. greater than one FB)
        figure(2);
        for i=1:num_of_episodes;
            line(rec_files.feeding_times_scored_hours_SyncedToEvtafTime(i,:), [0.98 0.98]);
        end
        
    end
    
    
end



%% OLD VERSION I TRIED TO DO WITH FOPEN BUT FAILED.
%
%
% for i=1:length(rec_files)
%     frecfile=fopen(rec_files(i).name,'r');
%     FB_found{i}=findstr(frecfile,'FB')
% end
%
%
%     fn{i}=rec_files(i).name;
%     ln=fgetl(fn(i));
%
%
% end


%% save everything

now_date=datestr(now,'ddmmmyyyy');
now_time=datestr(now,'HHMM');

mkdir(['lt_rec_files_find_FB_results_' now_date '_' now_time '_using' batch]);
cd(['lt_rec_files_find_FB_results_' now_date '_' now_time '_using' batch]);

saveas(figure(1), 'FB_times_possibly_with_feeding_time','fig')

% saveas(figure(2), 'song_times_possibly_with_feeding_time','fig')

save('rec_files_data','-struct','rec_files');

% warn to check using evsonganaly that sorting by FB=1 or FB>1
% distinguishes man trigger from real song

% disp('make sure real song has >1 FB and man trigger has only one -manually inspect these two batch files using evsonganaly');
%
% open(['../' batch '.rec_FB_MoreThanOne'])
% open(['../' batch '.rec_FB_OnlyOne']);
cd(directory);


% disp('save 2 new batch files which are RealSongsChecked and RealVocalizationChecked (any vocalization trigger, not just song) and rerun this script on those batch files')



