nhit=[0,0];nnt=[0,0];
for ii = 1:length(trigs)
	ind = 1;
	if (trigs(ii).ftime>=12)
		ind=2;
	end
	nhit(ind) = nhit(ind) + length(trigs(ii).ttimes);
	nnt(ind)  = nnt(ind)  + trigs(ii).totaltargetnt + length(trigs(ii).nontrignt);
end
disp(num2str(nhit));
disp(num2str(nnt));
disp(num2str(nhit./nnt));
