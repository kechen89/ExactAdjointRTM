BINDIR = ./../Bin
SRCDIR = ./../Src
FATPAR = labelfat=5 titlefat=5 plotfat=5

default: model_1d.ps frd1d.ps rtm1d.ps rt1dinv.ps rt1dinv_comp.ps rt1dinv_res.ps \
	 model_2d.ps rtm2d.ps

%.ps: %.v
	pstexpen $*.v $*.ps

${BINDIR}/%.x:
	! pushd ${SRCDIR} ; gmake ${BINDIR}/$*.x; popd
#
# 1-D Reverse-time Migration Example: 
#
#	wavelet.H 	: source wavelet to make initial exploding reflector wavefield
#	vel_1d.H    	: velocity model
#	init_1d.H 	: initial exploding reflector wavefield
#	data_1d.H    	: data with forward modeling
#	model_1d.H    	: image with exact adjoint RTM
#	model_1d_conv.H  : image with conventional RTM
#	rt1dinv.H    	: image with inversion
#
wavelet.H:
	Wavelet n1=1024 d1=0.001 wavelet=ricker2 phase=minimum tdelay=.00 | Window n1=512 > wavelet.H out=./wavelet.H@

vel_1d.H:
	Vel > v1.H vc=2500. n1=100 n2=1 n3=1 o1=0. d1=5. d2=5. d3=5.
	Vel > v2.H vc=4000. n1=100 n2=1 n3=1 o1=0. d1=5. d2=5. d3=5.
	Vel > v3.H vc=4750. n1=312 n2=1 n3=1 o1=0. d1=5. d2=5. d3=5.
	Merge space=n axis=1 v1.H v2.H v3.H > vel_1d.H out=./vel_1d.H@
	Rm v1.H v2.H v3.H

vel_1dn.H:
	Vel > v1.H vc=2750. n1=100 n2=1 n3=1 o1=0. d1=5. d2=5. d3=5.
	Vel > v2.H vc=4400. n1=100 n2=1 n3=1 o1=0. d1=5. d2=5. d3=5.
	Vel > v3.H vc=4750. n1=312 n2=1 n3=1 o1=0. d1=5. d2=5. d3=5.
	Merge space=n axis=1 v1.H v2.H v3.H > vel_1dn.H out=./vel_1dn.H@
	Rm v1.H v2.H v3.H

init_1d.H: wavelet.H
	Spike n1=512 n2=1 nsp=2 k1=100,200 mag=1.,0.7 > junk.H
	Filter < junk.H filter=wavelet.H pad=zero > init_1d.H out=./init_1d.H@
	echo " d1=5." >> init_1d.H
	Rm junk.H

data_1d.H: ${BINDIR}/test_rt1d.x init_1d.H vel_1d.H
	${BINDIR}/test_rt1d.x adj=0 snap=1 dt=.001 dz=5. nt=870 vel=vel_1d.H snapout=snap1d_frd.H < init_1d.H > data_1d.H out=./data_1d.H@
	
#--------------------------------
data_1dn.H: ${BINDIR}/test_rt1d.x data_1d.H
	Noise < data_1d.H type=1 var=0.03 > data_1dn.H
#--------------------------------
	
model_1d.H: ${BINDIR}/test_rt1d.x data_1d.H vel_1d.H
	${BINDIR}/test_rt1d.x adj=1 snap=1 dt=.001 dz=5. vel=vel_1d.H snapout=snap1d_rtm.H < data_1d.H > model_1d.H out=./model_1d.H@
	
model_1d_conv.H: ${BINDIR}/test_rt1d_conv.x data_1d.H vel_1d.H
	${BINDIR}/test_rt1d_conv.x adj=1 snap=0 dt=.001 dz=5. vel=vel_1d.H < data_1d.H > model_1d_conv.H out=model_1d_conv.H@
	
rt1dinv.H: ${BINDIR}/Rt1dinv.x data_1d.H vel_1d.H
	${BINDIR}/Rt1dinv.x niter=10 vel=vel_1d.H res=res_1d.H < data_1d.H > rt1dinv.H out=./rt1dinv.H@

#--------------------------------------------
rt1dinv_comp1.H: ${BINDIR}/Rt1dinv.x data_1dn.H vel_1d.H
	${BINDIR}/Rt1dinv.x niter=20 vel=vel_1d.H res=res_1d.H < data_1dn.H > rt1dinv_comp1.H out=./rt1dinv_comp1.H@

rt1dinv_comp2.H: ${BINDIR}/Rt1dinv.x data_1dn.H vel_1d.H
	${BINDIR}/Rt1dinv.x niter=20 vel=vel_1d.H mode=1 res=res_1d2.H < data_1dn.H > rt1dinv_comp2.H out=./rt1dinv_comp2.H@
#--------------------------------------------

#
# 2-D Reverse-time Migration Example: 
#
#	wavelet.H 	: source wavelet to make initial exploding reflector wavefield
#	vel_2d.H    	: velocity model
#	init_2d.H 	: initial exploding reflector wavefield
#	data_2d.H    	: data with forward modeling
#	model_2d.H    	: image with exact adjoint RTM
#	model_2d_conv.H  : image with conventional RTM
#	rt2dinv.H    	: image with inversion
#

vel_2d.H:
	Vel > v1.H vc=2000. n1=100 n2=800 n3=1 o1=0. o2=0. d1=5. d2=5. d3=5.
	Vel > v2.H vc=2600. n1=100 n2=800 n3=1 o1=0. o2=0. d1=5. d2=5. d3=5.
	Vel > v3.H vc=3200. n1=312 n2=800 n3=1 o1=0. o2=0. d1=5. d2=5. d3=5.
	Merge space=n axis=1 v1.H v2.H v3.H > vel_2d.H out=./vel_2d.H@
	Rm v1.H v2.H v3.H

vel_2ds.H:
	Vel > v1.H vc=2000. n1=100 n2=400 n3=1 o1=0. o2=0. d1=5. d2=5. d3=5.
	Vel > v2.H vc=2600. n1=100 n2=400 n3=1 o1=0. o2=0. d1=5. d2=5. d3=5.
	Vel > v3.H vc=3200. n1=312 n2=400 n3=1 o1=0. o2=0. d1=5. d2=5. d3=5.
	Merge space=n axis=1 v1.H v2.H v3.H > vel_2ds.H out=./vel_2ds.H@
	Rm v1.H v2.H v3.H

init_2d.H: wavelet.H init_1d.H
	Spike k1=-1 n1=512 n2=20 n3=1 d1=5. d2=5. > junk.H
	Pad < init_1d.H n2=760 extend=1 > junk1.H
	Merge space=no axis=2 junk.H junk1.H junk.H > junk2.H
	Spike n1=512 n2=800 n3=1 nsp=7 k1=300,300,300,300,300,300,300  k2=397,398,399,400,401,402,403 mag=0.4,0.5,0.66,1.,0.66,0.5,0.4 > junk.H
	Filter < junk.H filter=wavelet.H pad=zero > junk1.H
	Add junk2.H junk1.H > init_2d.H out=./init_2d.H@
	Rm junk.H junk1.H junk2.H

init_2ds.H: wavelet.H init_1d.H
	Spike k1=-1 n1=512 n2=20 n3=1 d1=5. d2=5. > junk.H
	Pad < init_1d.H n2=760 extend=1 > junk1.H
	Merge space=no axis=2 junk.H junk1.H junk.H > junk2.H
	Spike n1=512 n2=800 n3=1 nsp=7 k1=300,300,300,300,300,300,300  k2=397,398,399,400,401,402,403 mag=0.4,0.5,0.66,1.,0.66,0.5,0.4 > junk.H
	Filter < junk.H filter=wavelet.H pad=zero > junk1.H
	Add junk2.H junk1.H > init_2d.H out=./init_2d.H@
	Rm junk.H junk1.H junk2.H

data_2d.H: ${BINDIR}/test_rt2d.x vel_2d.H init_2d.H
	echo "this will take about 3 minutes on a PC with intel Pentium-4 3.0Ghz CPU" 
	${BINDIR}/test_rt2d.x adj=0 dt=.001 nt=1250 vel=vel_2d.H < init_2d.H > data_2d.H out=./data_2d.H@

#--------------------------------
data_2dn.H: ${BINDIR}/test_rt2d.x data_2d.H
	Noise < data_2d.H type=1 var=0.03 > data_2dn.H
#--------------------------------
	
model_2d.H: ${BINDIR}/test_rt2d.x vel_2d.H data_2d.H
	echo "this will take about 4 minutes on a PC with intel Pentium-4 3.0Ghz CPU" 
	${BINDIR}/test_rt2d.x adj=1 vel=vel_2d.H < data_2d.H > model_2d.H out=./model_2d.H@

model_2d_inv.H: ${BINDIR}/test_rt2d.x vel_2d.H data_2d.H
	echo "this will take about 4 minutes on a PC with intel Pentium-4 3.0Ghz CPU" 
	${BINDIR}/test_rt2dinv.x adj=1 vel=vel_2d.H niter=10 res=res2d.H < data_2d.H > model_2d_inv.H

model_2d_conv.H: ${BINDIR}/test_rt2d_conv.x vel_2d.H data_2d.H 
	echo "this will take about 4 minutes on a PC with intel Pentium-4 3.0Ghz CPU" 
	${BINDIR}/test_rt2d_conv.x adj=1 vel=vel_2d.H < data_2d.H > model_2d_conv.H out=./model_2d_conv.H@

rt2dinv.H: ${BINDIR}/Rt2dinv.x vel_2d.H data_2d.H
	echo "this will take about 75 minutes hours on a PC with intel Pentium-4 3.0Ghz CPU" 
	date
	${BINDIR}/Rt2dinv.x vel=vel_2d.H niter=10 res=res_2d.H < data_2d.H > rt2dinv.H out=./rt2dinv.H@
	date

#--------------------------------------------------------------
rt2dinv_comp1.H: ${BINDIR}/Rt2dinv.x vel_2d.H data_2dn.H
	echo "this will take about 75 minutes on a PC with intel Pentium-4 3.0Ghz CPU" 
	date
	${BINDIR}/Rt2dinv.x vel=vel_2d.H niter=10 mode=0 res=res_2d_comp1.H < data_2dn.H > rt2dinv_comp1.H out=./rt2dinv_comp1.H@
	date

rt2dinv_comp2.H: ${BINDIR}/Rt2dinv.x vel_2d.H data_2dn.H
	echo "this will take about 75 minutes hours on a PC with intel Pentium-4 3.0Ghz CPU" 
	date
	${BINDIR}/Rt2dinv.x vel=vel_2d.H niter=10 mode=1 res=res_2d_comp2.H < data_2dn.H > rt2dinv_comp2.H out=./rt2dinv_comp2.H@
	date
#--------------------------------------------------------------
#
#   Making figures
#

#  Fig. 3 in the paper
model_1d.v: vel_1d.H init_1d.H
	sfgraph < init_1d.H title="Initial Wavefield" label1="Depth(m)" label2="Amplitude" ${FATPAR} labelsz=10 titlesz=10 parallel2=n > a.v
	sfgraph < vel_1d.H title="Velocity Model" label1="Depth(m)" label2="velocity(m/s)" ${FATPAR} labelsz=10 titlesz=10 min2=2000. max2=5000. parallel2=n >b.v
	vp_OverUnderAniso a.v b.v > model_1d.v
	/bin/rm a.v b.v

#  Fig. 4 in the paper
frd1d.v: data_1d.H
	Dots < data_1d.H wantaxis=0 ${FATPAR} labels=" " title=" "> /dev/null out=a.v
	vppen a.v align=lb > b.v
	vppen xsize=4.35 ysize=1. < b.v > a.v
	vppen xshift=1.15 yshift=5.0 < a.v > b.v
	echo "o1=0.0 d1=5.0 o2=0.0 d2=0.001 n2=868" >> snap1d_frd.H
	sfgrey < snap1d_frd.H pclip=100 title=" " ${FATPAR} label1="Depth(m)" wherexlabel=b label2="Time(sec)" labelsz=7 labelfat=3 parallel2=n > a.v
	vppen a.v align=lb > c.v
	vppen xsize=5. ysize=5. < c.v > a.v
	vppen xshift=.5 < a.v > c.v
	vppen erase=o b.v c.v > frd1d.v
	/bin/rm a.v b.v c.v

#  Fig. 5 in the paper
rtm1d.v: model_1d.H
	Dots < model_1d.H wantaxis=0 ${FATPAR} labels=" " transp=1 title=" "> /dev/null out=a.v
	vppen a.v align=lb > b.v
	vppen xsize=0.8 ysize=4.5 < b.v > a.v
	vppen xshift=.3 yshift=0.0 < a.v > b.v
	Reverse < snap1d_rtm.H > junk.H
	Window < junk.H f2=0 j2=2 > junk1.H
	echo "o1=0.0 d1=5.0 o2=0.0 d2=0.002" >> junk1.H
	sfgrey < junk1.H pclip=99 title=" " ${FATPAR} label1="Depth(m)" label2="Time(sec)" labelsz=7 labelfat=3 whereylabel=r grid1=0 grid2=0  labelrot=y parallel2=n > a.v
	vppen a.v align=lb > c.v
	vppen xsize=5. ysize=5. < c.v > a.v
	vppen xshift=1.4 < a.v > c.v
	vppen erase=o b.v c.v > rtm1d.v
	/bin/rm a.v b.v c.v 
	Rm junk.H junk1.H

#  Fig. 6 in the paper
rt1dinv.v: init_1d.H model_1d.H model_1d_conv.H rt1dinv.H
	Merge space=0 axis=2 init_1d.H model_1d.H model_1d_conv.H rt1dinv.H >junk.H
	Dots < junk.H yreverse=1 title=" " ${FATPAR} labels="(a) INPUT:(b) RTM(Adjoint):(c) RTM(Conventional):(d) LS RTM" labelsz=6 label1="Depth(m)" plotfat=4 > /dev/null out=rt1dinv.v
	Rm junk.H

#  Fig. 7 in the paper
rt1dinv_comp.v: data_1dn.H rt1dinv_comp1.H rt1dinv_comp2.H
	Dots < data_1dn.H yreverse=1 title=" " ${FATPAR} labels="(a) Noisy input" labelsz=20 label1="Time (sec)" > /dev/null out=a.v
	vppen xsize=4.2 ysize=3. < a.v > b.v
	vppen yshift=0.2 < b.v > a.v
	Dots < rt1dinv_comp1.H yreverse=1 title=" " ${FATPAR} labels="(b) Exact adjoint pair" labelsz=20 label1="Depth (m)" > /dev/null out=b.v
	vppen xsize=4.2 ysize=3. < b.v > c.v
	vppen yshift=0.2 < c.v > b.v
	Dots < rt1dinv_comp2.H yreverse=1 title=" " ${FATPAR} labels="(c) Conventional adjoint \n    pair" labelsz=20 label1="Depth (m)" > /dev/null out=c.v
	vppen xsize=4.2 ysize=3. < c.v > d.v
	vppen yshift=0.2 < d.v > c.v
	vp_OverUnderAniso a.v b.v c.v > rt1dinv_comp.v
	/bin/rm a.v b.v c.v

#  Fig. 8 in the paper
rt1dinv_res.v: rt1dinv_comp1.H rt1dinv_comp2.H
	Merge space=0 axis=2 res_1d.H res_1d2.H >junk.H
	sfgraph < junk.H labelsz=6 title=" " ${FATPAR} label1="Iteration" label2="Mean square errors" dash=0,1 legend=y curvelabel="Exact adjoint pair LSRTM: Conventional adjoint pair LSRTM" parallel2=n >  junk.v
	vp_annotate < junk.v batch=y text=text_file > rt1dinv_res.v
	Rm junk.H
	rm junk.v

#  Fig. 9 in the paper
model_2d.v: data_2d.H vel_2d.H
	sfgrey < init_2d.H titlesz=11 labelsz=11 ${FATPAR} label1="Depth(m)" label2="Surface(m)" title="(a) Source Wavefield" parallel2=n >a.v
	sfgrey < vel_2d.H titlesz=11 labelsz=11 ${FATPAR} pclip=98 bias=1600. allpos=y title="(b) Velocity Model" label1="Depth(m)" label2="Surface(m)" parallel2=n > b.v
	sfgrey < data_2d.H titlesz=11 labelsz=11 ${FATPAR} title="(c) Modeled data" label1="Time(s)" label2="Surface(m)" parallel2=n > c.v
	vp_SideBySideAniso a.v b.v c.v > model_2d.v
	/bin/rm a.v b.v c.v

#  Fig. 10 in the paper
rtm2d.v: model_2d.H model_2d_conv.H rt2dinv.H
	sfgrey < model_2d.H titlesz=11 ${FATPAR} labelsz=11 title="(a) RTM (Adjoint)" label1="Depth(m)" label2="Surface(m)" parallel2=n > a.v
	sfgrey < model_2d_conv.H titlesz=11 ${FATPAR} labelsz=11 title="(b) RTM (Conventional)" label1="Depth(m)" label2="Surface(m)" parallel2=n > b.v
	sfgrey < rt2dinv.H titlesz=11 ${FATPAR} labelsz=11 title="(c) LS RTM" label1="Depth(m)" label2="Surface(m)" parallel2=n > c.v
	vp_SideBySideAniso a.v b.v c.v > rtm2d.v
	/bin/rm a.v b.v c.v

clean:
	Rm *.H
	/bin/rm a.out core *.v *.ps junk* xyz*
