
%designed to call pitchsyntaxanal, and then inactivanal.m
%extra baseline day???

clear date
clear muanal
clear avls
clear ac
clear graphvals pathvl catchvl usex  wn wnrev
sumpath='/oriole4/bk20bk45/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
avls.initmean{1}=3319
avls.initsd{1}=58.5;
avls.maxmeans=[4731 4765]

avls.initmean{2}=5145
avls.initsd{2}=95.24


avls.initmean{3}=5145
avls.initsd{3}=95.24

avls.initmean{4}=3237
avls.initsd{4}=50.1
extractvals=1
avls.changenote=0

%index into cell array is the note
graphvals.wnon=1;
wn(1).tmon{1}={'2008-7-13 7' '2008-8-10 7'}
wn(1).tmoff{1}={'2008-8-07 7' '2008-8-12 7'}
%turned off 7-24


wn(4).tmon=wn(1).tmon;
wn(4).tmoff=wn(1).tmoff;

wn(2).tmon{1}={'2008-7-13 7' '2008-7-17 7'} 
wn(2).tmon{2}={'2008-7-26 7' '2008-7-28 7'}
wn(2).tmoff{1}={'2008-7-17 7' '2008-7-22 7'}
 wn(2).tmoff{2}=   {'2008-7-28 7' '2008-8-3 7' }

wn(3).tmon=wn(2).tmon;
wn(3).tmoff=wn(2).tmoff;

% wnrev(3).tmon{1}={'2008-7-23 7'}
% wnrev(3).tmoff{1}={'2008-7-26 7'}

wnrev(3).tmon{1}={'2008-7-23 7'}
wnrev(3).tmoff{1}={'2008-7-26 7'}

wnrev(3).tmon{2}={'2008-8-03 7'}
wnrev(3).tmoff{2}={'2008-8-07 7'}

wnbnds{1}=[3000 4000]
wnbnds{2}=[4000 6000]

edges{1}=[3000:100:4800]
edges{2}=[4600:75:5800]
edges{3}=[4600:75:5800]
edges{4}=edges{1}
edges{5}=edges{1}

% 6/5, 6/6 and 6/7 again 6/14. 6/15, 6/16
% wn{1}.freqlo=[3000 3000 3000 3000 3000     3745 3600 3470 3245]
% wn{1}.freqhi=[3400 3450 3500 3600 3675     4000 4000 4000 4000 ]
% wn{2}.freqlo=[5100 5000  4950 4940 4860     4000 4000 4000 4000]
% wn{2}.freqhi=[6000 6000  6000 6000 6000     4950 4900 4875 5220]
% 
% wn{1}.tmon={'2008-06-05 07:00' '2008-06-05 13:15' '2008-06-05 16:17' '2008-06-06 10:20' '2008-06-07 07:00:00' '2008-06-14 07:00' '2008-06-15 07:00' '2008-06-16 07:00' '2008-06-21 07:00'};
% wn{1}.tmoff=[wn{1}.tmon(2:4) '2008-06-14 07:00' wn{1}.tmon(6:8) '2008-06-16 21:00' '2008-07-02 21:00']
% wn{2}.tmon=wn{1}.tmon;
% wn{2}.tmoff=wn{1}.tmoff;
% 
% wn{3}=wn{2}
% wn{4}=wn{1}
% wn{5}=wn{1};

clear pathvl vals tmn tmf




pathvl{1}='/oriole4/bk20bk45/probein/'
catchvl{1}='batch10.keep.rand'
tmn{1}='07:00'
tmf{1}='09:15'

% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole4/bk20bk45/lid710/'
catchvl{2}='batchcomb'
tmn{2}='09:15'
tmf{2}='12:30'

pathvl{3}='/oriole4/bk20bk45/ac710/'
catchvl{3}='batch.keep.rand'
tmn{3}='12:30'
tmf{3}='16:00'

pathvl{4}='/oriole4/bk20bk45/750mu710/'
catchvl{4}='batch.keep'
tmn{4}='16:00'
tmf{4}='21:00'

pathvl{5}='/oriole4/bk20bk45/ac710-2/'
catchvl{5}='batch11.keep'
tmn{5}='07:00'
tmf{5}='10:13'

pathvl{6}='/oriole4/bk20bk45/500mu711/'
catchvl{6}='batch'
tmn{6}='10:13'
tmf{6}='14:13'

pathvl{7}='/oriole4/bk20bk45/ac711/'
catchvl{7}='batch11.keep.rand'
tmn{7}='14:13'
tmf{7}='21:00'

pathvl{8}='/oriole4/bk20bk45/ac721/'
catchvl{8}='batch.catch.keep'
tmn{8}='14:20'
tmf{8}='21'

pathvl{9}='/oriole4/bk20bk45/ac723/'
catchvl{9}='batch23.keep.catch'
tmn{9}='16:50'
tmf{9}='21'

pathvl{10}='/oriole4/bk20bk45/ac712/'
catchvl{10}='batch.keep.rand'
tmn{10}='11:34'
tmf{10}='14:42'

pathvl{11}='/oriole4/bk20bk45/600mu712/'
catchvl{11}='batch'
tmn{11}='14:42'
tmf{11}='18:12'


pathvl{12}='/oriole4/bk20bk45/wnon713/'
catchvl{12}='batchcomb'
tmn{12}='7'
tmf{12}='12:30'

pathvl{13}='/oriole4/bk20bk45/600715/'
catchvl{13}='batch.catch'
tmn{13}='12:07'
tmf{13}='16:07'

pathvl{14}='/oriole4/bk20bk45/ac715/'
catchvl{14}='batch16'
tmn{14}='7'
tmf{14}='21'

pathvl{15}='/oriole4/bk20bk45/ac715/'
catchvl{15}='batch17.catch.keep'
tmn{15}='7'
tmf{15}='11:07'

pathvl{16}='/oriole4/bk20bk45/600717/'
catchvl{16}='batch.catch.keep'
tmn{16}='11:07'
tmf{16}='14:34'

pathvl{17}='/oriole4/bk20bk45/ac717/'
catchvl{17}='batch.catch.keep'
tmn{17}='14:34'
tmf{17}='21'

pathvl{18}='/oriole4/bk20bk45/ac715/'
catchvl{18}='batch15.catch.keep'
tmn{18}='16:07'
tmf{18}='21'

pathvl{19}='/oriole4/bk20bk45/ac712-2/'
catchvl{19}='batch12.keep'
tmn{19}='18:12'
tmf{19}='21'


pathvl{20}='/oriole4/bk20bk45/ac717/'
catchvl{20}='batch18.catch.keep'
tmn{20}='7'
tmf{20}='9:13'


pathvl{21}='/oriole4/bk20bk45/700718/'
catchvl{21}='batch18.catch.keep'
tmn{21}='9:13'
tmf{21}='13:15'

pathvl{22}='/oriole4/bk20bk45/ac718/'
catchvl{22}='batch.catch.keep'
tmn{22}='12:55'
tmf{22}='21:00'


pathvl{23}='/oriole4/bk20bk45/ac718/'
catchvl{23}='batch19.catch.keep'
tmn{23}='7'
tmf{23}='11:28'

%actually the 19th.
pathvl{24}='/oriole4/bk20bk45/800718/'
catchvl{24}='batch.catch.keep'
tmn{24}='11:28'
tmf{24}='15:15'


pathvl{25}='/oriole4/bk20bk45/canin/'
catchvl{25}='batch.keep'
tmn{25}='7'
tmf{25}='21:00'

pathvl{26}='/oriole4/bk20bk45/ac719/'
catchvl{26}='batch21.catch.keep'
tmn{26}='7'
tmf{26}='10:20'

pathvl{27}='/oriole4/bk20bk45/800721/'
catchvl{27}='batch21.catch.keep'
tmn{27}='10:20'
tmf{27}='14:20'

pathvl{27}='/oriole4/bk20bk45/800721/'
catchvl{27}='batch21.catch.keep'
tmn{27}='10:20'
tmf{27}='14:20'

pathvl{28}='/oriole4/bk20bk45/wnrev721/'
catchvl{28}='batch.keep.catch'
tmn{28}='7'
tmf{28}='12:45'

pathvl{29}='/oriole4/bk20bk45/900723/'
catchvl{29}='batch.keep'
tmn{29}='12:45'
tmf{29}='16:50'

pathvl{30}='/oriole4/bk20bk45/ac723/'
catchvl{30}='batchcomb24'
tmn{30}='7'
tmf{30}='11:00'

pathvl{31}='/oriole4/bk20bk45/900724/'
catchvl{31}='batchcomb'
tmn{31}='11:00'
tmf{31}='15:00'

pathvl{32}='/oriole4/bk20bk45/ac724/'
catchvl{32}='batch24.keep.rand'
tmn{32}='15:00'
tmf{32}='21'

% pathvl{33}='/oriole4/bk20bk45/ac724/'
% catchvl{33}='batch25.keep.rand'
% tmn{33}='7'
% tmf{33}='21'

% pathvl{34}='/oriole4/bk20bk45/wnon2-b/'
% catchvl{34}='batch27.catch.keep'
% tmn{34}='7'
% tmf{34}='21'

pathvl{35}='/oriole4/bk20bk45/wnon2-b/'
catchvl{35}='batch.catch.keep'
tmn{35}='7'
tmf{35}='9:53'


pathvl{36}='/oriole4/bk20bk45/875_728/'
catchvl{36}='batch.keep.catch'
tmn{36}='9:53'
tmf{36}='13:56'

pathvl{37}='/oriole4/bk20bk45/ac730/'
catchvl{37}='batch01.catch.keep'
tmn{37}='7'
tmf{37}='14:00'



pathvl{38}='/oriole4/bk20bk45/lido730/'
catchvl{38}='batchcomb'
tmn{38}='12:04'
tmf{38}='16:17'

pathvl{39}='/oriole4/bk20bk45/ac728/'
catchvl{39}='batch28.catch.keep'
tmn{39}='13:56'
tmf{39}='21'


pathvl{40}='/oriole4/bk20bk45/ac801/'
catchvl{40}='batch02.catch.keep'
tmn{40}='7'
tmf{40}='09:45'

pathvl{41}='/oriole4/bk20bk45/875_802/'
catchvl{41}='batch.catch.keep'
tmn{41}='09:45'
tmf{41}='13:20'

pathvl{42}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{42}='batch03mu.catch.keep'
tmn{42}='11:48'
tmf{42}='16:02'

pathvl{43}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{43}='batch03ac.catch.keep'
tmn{43}='16:02'
tmf{43}='21'

pathvl{44}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{44}='batch04ac1.catch.keep'
tmn{44}='7'
tmf{44}='11:24'

pathvl{45}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{45}='batchcomb04mu'
tmn{45}='11:24'
tmf{45}='15:04'

pathvl{46}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{46}='batch04ac2.catch.keep'
tmn{46}='15:04'
tmf{46}='21'

pathvl{47}='/oriole4/bk20bk45/wnswitch2/'
catchvl{47}='batch03.keep.catch'
tmn{47}='7'
tmf{47}='11:48'


pathvl{48}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{48}='batch05ac1.catch.keep'
tmn{48}='7'
tmf{48}='11:41'

pathvl{49}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{49}='batchcomb05mu'
tmn{49}='11:41'
tmf{49}='15:12'

pathvl{50}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{50}='batch05ac2.catch.keep'
tmn{50}='15:12'
tmf{50}='21'

pathvl{51}='/oriole4/bk20bk45/wnswitch2lidon/'
catchvl{51}='batch06.catch.keep'
tmn{51}='7'
tmf{51}='21'

pathvl{52}='/oriole4/bk20bk45/ac802/'
catchvl{52}='batch.catch.keep'
tmn{52}='13:20'
tmf{52}='21'


pathvl{53}='/oriole4/bk20bk45/ac719/'
catchvl{53}='batch19.catch.keep'
tmn{53}='15:15'
tmf{53}='21'

pathvl{54}='/oriole4/bk20bk45/ac728/'
catchvl{54}='batch30.catch.keep'
tmn{54}='7'
tmf{54}='12:04'

pathvl{55}='/oriole4/bk20bk45/ac730/'
catchvl{55}='batch30.keep'
tmn{55}='17'
tmf{55}='21'

pathvl{56}='/oriole4/bk20bk45/806_postdata/'
catchvl{56}='batch26.keep.rand'
tmn{56}='7'
tmf{56}='21'

% pathvl{57}='/oriole4/bk20bk45/wnon/'
% catchvl{57}='batch.catch.keep'
% tmn{57}='7'
% tmf{57}='21'

pathvl{58}='/oriole4/bk20bk45/wnon_newtemp/'
catchvl{58}='batch.catch.keep'
tmn{58}='7'
tmf{58}='21'

pathvl{59}='/oriole4/bk20bk45/wnon_newtemp/'
catchvl{59}='batch14.catch.keep'
tmn{59}='7'
tmf{59}='21'

pathvl{60}='/oriole4/bk20bk45/wnon_newtemp/'
catchvl{60}='batch15.catch.keep'
tmn{60}='7'
tmf{60}='21'

pathvl{61}='/oriole4/bk20bk45/wnon_newtemp2/'
catchvl{61}='batch16.catch.keep'
tmn{61}='7'
tmf{61}='21'


pathvl{62}='/oriole4/bk20bk45/wnon_newtemp3/'
catchvl{62}='batch19.keep.rand'
tmn{62}='7'
tmf{62}='21'

pathvl{63}='/oriole4/bk20bk45/wnoff919/'
catchvl{63}='batch23.keep.rand'
tmn{63}='7'
tmf{63}='21'

pathvl{64}='/oriole4/bk20bk45/wnon1024/'
catchvl{64}='batch26.keep.catch'
tmn{64}='7'
tmf{64}='21'

pathvl{65}='/oriole4/bk20bk45/wnoff919/'
catchvl{65}='batch25.keep.rand'
tmn{65}='7'
tmf{65}='21'

pathvl{66}='/oriole4/bk20bk45/wn1121/'
catchvl{66}='batchcomb'
tmn{66}='7'
tmf{66}='21'

pathvl{67}='/oriole4/bk20bk45/tmptest/'
catchvl{67}='batch.keep.rand'
tmn{67}='7'
tmf{67}='21'

pathvl{68}='/oriole4/bk20bk45/wnoff928/'
catchvl{68}='batch28.keep'
tmn{68}='7'
tmf{68}='21'

pathvl{69}='/oriole4/bk20bk45/wnoff928/'
catchvl{69}='batch01.keep'
tmn{69}='7'
tmf{69}='21'

pathvl{70}='/oriole4/bk20bk45/wnon713/'
catchvl{70}='batch14.catch'
tmn{70}='7'
tmf{70}='21'

pathvl{71}='/oriole4/bk20bk45/ac719/'
catchvl{71}='batch20.catch.keep'
tmn{71}='7'
tmf{71}='21'

pathvl{72}='/oriole4/bk20bk45/wnrev721/'
catchvl{72}='batch22.catch'
tmn{72}='7'
tmf{72}='21'

pathvl{73}='/oriole4/bk20bk45/wnon2-b/'
catchvl{73}='batch27.catch'
tmn{73}='7'
tmf{73}='21'

pathvl{74}='/oriole4/bk20bk45/ac728/'
catchvl{74}='batch29.catch'
tmn{74}='7'
tmf{74}='21'



usex(1:39)=0;
ac(39)=[0];
avls.analind=[1:39]


date{39}='';

numnotes=2
notes='ab'

fbins={};
tbinshft{1}=0.015;
  tbinshft{2}=0.012;
tbinshft{3}=0.01;
  tbinshft{4}=0.01;
  tbinshft{5}=0.01
NFFT(1)=256
NFFT(2)=256
 NFFT(3)=256
NFFT(4)=512
 NFFT(5)=512

fbins{1}=[3000,4500];
fbins{2}=[4600,5800];
fbins{3}=[4600,5800];
fbins{4}=[3000,4500];
fbins{5}=[3000,4500];

clear NT    
NT{1}='a'
 NT{2}='b'
NT{3}='b'
NT{4}='a'
% NT{5}='a'
% NT{3}='g'

PRENT{1}='--';PSTNT{1}='';
PRENT{2}='b';PSTNT{2}='';
PRENT{3}='';PSTNT{3}='b';
PRENT{4}='a';PSTNT{4}='a';

avls.exclude=[0 1 0 1]

% PRENT{5}='b';PSTNT{5}='';
  %%%plotvals
acon=[5 7 8 9 10 12 14 15 17 18 19 20 22 23 25 26 28 30 32  35 37 39 40 43 44 46 47 48 50 52:55 70:74] ;  
muon=[6 11 13 16 21 24 27 29 31 36 38 41 42 45 49 ]
avls.muanal{3}{1}=[6 11 13 16 21 24 27 ];
avls.muanal{3}{2}=[36 38 41];
avls.revanal{3}{1}=[29 31];
avls.revanal{3}{2}=[42 45 49];
avls.revanal{1}=[];

avls.muanal{1}{1}=[6 11 13 16 21 24 ];


diron=[] 
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon

avls.edges=edges
graphvals.numcol=2
graphvals.col='rk'
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.tickht=8000
graphvals.edges{1}=[6700:80:7900]
graphvals.edges{2}=[2000:40:2500]
graphvals.edges{3}=[1300:40:1800]
graphvals.colvals=colvals
graphvals.pltpnts=1

% graphvals.timon=timon;
% graphvals.timoff=timoff;
graphvals.acon=ac;
graphvals.date=date;

graphvals.wn=wn;
graphvals.wnrev=wnrev;
avls.diron=diron;
graphvals.colvals=colvals
graphvals.chunkdata=1

avls.pvls=pathvl
avls.datfile=matfilename;
avls.cvl=catchvl
avls.sumpath=sumpath
avls.mtflnm=matfilename  
avls.supanal=0
avls.NT=NT
avls.NFFT=NFFT
avls.fbins=fbins
avls.tshft=tbinshft
avls.usex=usex
avls.tmon=tmn;% NT{4}='a'
% NT{5}='a'

avls.tmoff=tmf;
avls.date=date;
avls.numnotes=numnotes
avls.bname='bk20bk45'
avls.mkfv=[35]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.muoffset=1.5/24;
avls.offset=1;
avls.acoffset=1.5/24;
avls.deadtm=.0833/24;
avls.acon=acon;
avls.muon=muon;
avls.pretm=3/24

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


