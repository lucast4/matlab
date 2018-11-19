% SaveDir='/bluejay4/lucas/Dropbox/SCIENCE/BRAINARD_LAB/PROJECTS/MATLAB';
% SaveDir = '/bluejay5/egret_data/lucas/Transferfiles';
SaveDir = '/home/lucast4/Dropbox/SCIENCE/BRAINARD_LAB/FIGURES/MATLAB';
timestamp=lt_get_timestamp(0);
randnum=randi(10000, 1);
SavePath=[SaveDir '/' timestamp '_rand' num2str(randnum) '.eps'];


% print('-depsc2', '-noui', '-adobecset', '-painters', SavePath)
% good:
% print('-depsc2', '-noui', '-painters', SavePath)
% export_fig(SavePath, '-pdf', '-transparent', '-depsc');  doesn't work ends up as bmp
print('-dpdf', '-painters', SavePath);
%% Using export_fig, poptentialyl better?

% -transparent

if (0) % not working so far. might need to download ghostscript
savetemp=[SaveDir '/export_fig_version.eps'];


export_fig('-painters', '-depsc', '-eps', '-CMYK', savetemp);
export_fig savetemp -png
end