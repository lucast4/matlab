
ps.minx=3
ps.maxx=4
ps.ac_col='k'
ps.mu_col=[0.4 0.4 1]

ps.flip=1
ps.plotraw=1
ps.addzero=0
ps.plotsep=1;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.plotpct=1
figure
ps.subdayinterp=1

[outvlaczcomb,outvlmuzcomb]=plotrevdynamics(sumdynrev(1:6),ps)
