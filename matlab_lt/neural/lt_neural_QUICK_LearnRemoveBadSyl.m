function sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, swnum, syltoken)
%% tells you whether you should keep a given syl for a given switch, learing
% 
%% syllabels that follow the target (on same motif up to 2+ syls) should be excluded

SylsBad = {...,
    {'pu69wh78', 'RALMANOvernightLearn1', 1, 'aabh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 6, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 6, 'jjbh(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 7, 'jjb(h)'}, ...
    {'pu69wh78', 'RALMANOvernightLearn1', 7, 'jjbh(h)'}, ...
    {'wh44wh39', 'RALMANlearn1', 2, 'cb(b)'}, ...
    {'wh44wh39', 'RALMANlearn3', 1, 'nh(h)'}, ...
    {'wh44wh39', 'RALMANlearn4', 2, 'dkc(c)'}, ...
};


%% ========== ask whether the input matches any of the bad syls
sylbad = 0;
for j=1:length(SylsBad)
   tmp1 = strcmp(SylsBad{j}{1}, bname);
   tmp2 = strcmp(SylsBad{j}{2}, ename);
   tmp3 = SylsBad{j}{3}==swnum;
   tmp4 = strcmp(SylsBad{j}{4}, syltoken);
   
   if all([tmp1 tmp2 tmp3 tmp4])
       sylbad=1;
       break
   end
end

