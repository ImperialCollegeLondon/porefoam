# converting porefoam two-phase flow results (grads.csv) to 
# relative permeability and capillary pressures 
# developed by Ali Q Raeini



def groupTime(xxx,indx,nGroup):
	xx=[row[indx] for row in xxx]
	xx_a1=np.zeros(int(len(xx)/nGroup))
	n_a1=nGroup - 2
	if n_a1>1:
		for i in range(0,len(xx_a1)):
			ixx2= (i+1)*nGroup
			group=xx[ixx2 - nGroup :ixx2]
			xx_a1[i]=(sum(group) - max(group) - min(group)) / n_a1 #filter out min and maxs
	else:
		for i in range(len(xx_a1)):
			ixx2= (i+1)*nGroup
			group=xx[ixx2 - nGroup :ixx2]
			xx_a1[i]=np.mean(group)
	return xx_a1



def groupTimeAvg(xxx=None,Names=None,VarName=None,inds=None,delNGroup=None):
	weight=groupTime(xxx,Names.index('S'+str(inds[0])+'-vol'),delNGroup)*1e18
	xx_a1SumWeight=weight
	xx_a1Sum=(weight*groupTime(xxx,Names.index('S'+str(inds[0])+'-'+VarName),delNGroup))
	for i in range(1,len(inds)):
		weight=groupTime(xxx,Names.index('S'+str(inds[i])+'-vol'),delNGroup)*1e18
		xx_a1SumWeight=xx_a1SumWeight + weight
		xx_a1Sum=xx_a1Sum + (weight*groupTime(xxx,Names.index('S'+str(inds[i])+'-'+VarName),delNGroup))
	xx_a1Avg=xx_a1Sum / xx_a1SumWeight
	return xx_a1Avg



# def groupTimeMin(xxx=None,Names=None,VarName=None,inds=None,delNGroup=None):
	# indx=Names.index('S'+str(inds[0])+'-'+VarName)
	# xx_aMin=min(groupTime(xxx,indx,delNGroup))
	# for i in range(1,len(inds)):
		# indx=Names.index('S'+str(inds[i])+'-'+VarName)
		# xx_aMin=min(xx_aMin,min(groupTime(xxx,indx,delNGroup).tolist()))
	# return xx_aMin



# def groupTimeMax(xxx=None,Names=None,VarName=None,inds=None,delNGroup=None):
	# indx=Names.index('S'+str(inds[0])+'-'+VarName)
	# xx_a1Max=max(groupTime(xxx,indx,delNGroup))
	# for i in range(1,len(inds)):
		# indx=Names.index(('S'+str(inds[i])+'-'+VarName))
		# xx_a1Max=max(xx_a1Max,max(groupTime(xxx,indx,delNGroup)))
	# return xx_a1Max



# def groupTimeSum(xxx=None,Names='',VarName='',inds=None,delNGroup=None):
	# indx=Names.index('S'+str(inds[0])+'-'+VarName)
	# xx_a1Sum=groupTime(xxx,indx,delNGroup)
	# for i in range(1,len(inds)):
		# indx=Names.index(('S'+str(inds[i])+'-'+VarName))
		# xx_a1Sum=(xx_a1Sum + groupTime(xxx,indx,delNGroup))
	# return xx_a1Sum






def gradsToRelPermSwPcKrwKroT(Names=[''],volInt=[[]],SinglePhaseRes=[[]],inds=[],nGroup=10):
	'''relative conductivities assuming same viscsity for single-phase flow and two-phase flow simulations'''
	#nGroup=10
	#Names = genfromtxt('grads_hdr.csv', delimiter=' ',dtype='string').tolist()
	#volInt = genfromtxt('grads.csv', delimiter=',',dtype='double',skip_header=1)[1:-1]
	#SinglePhaseRes=[[]]
	#SinglePhaseRes[0] = genfromtxt('grads_SP.csv', delimiter=',',dtype='double')

	# compute integrals
	t=groupTime(volInt,0,nGroup)
	alpha1=groupTimeAvg(volInt,Names,('alpha'),inds,nGroup)
	Pc=groupTimeAvg(volInt,Names,('Pc'),inds,nGroup)
	UzAvg=groupTimeAvg(volInt,Names,('U'),inds,nGroup)
	f_1=groupTimeAvg(volInt,Names,('f_1'),inds,nGroup)
	Uz1Avg=(UzAvg*f_1)
	Uz2Avg=UzAvg - Uz1Avg
	viscEV=groupTimeAvg(volInt,Names,('viscE'),inds,nGroup)
	viscEV1=groupTimeAvg(volInt,Names,('viscE_1'),inds,nGroup)
	viscEV2=viscEV - viscEV1 - 1e-6
	UzAvgSP=groupTimeAvg(SinglePhaseRes,Names,'U',inds,1)
	viscEVSP=groupTimeAvg(SinglePhaseRes,Names,(('viscE')),inds,1)
	viscEVSP=viscEVSP[-1]
	SwPcKrs=[[],[],[],[],[]]
	SwPcKrs[0]=(alpha1)
	SwPcKrs[1]=(Pc)
	SwPcKrs[2]=((Uz1Avg*Uz1Avg) / (viscEV1-1e-16) / ((UzAvgSP*UzAvgSP) / viscEVSP))
	SwPcKrs[3]=((Uz2Avg*Uz2Avg) / (viscEV2-1e-16) / ((UzAvgSP*UzAvgSP) / viscEVSP))
	SwPcKrs[4]=(t)
	return SwPcKrs




def gradsToRelPermSwPcKrwKroTRev(Names=[''],volInt=[[]],SinglePhaseRes=[[]],inds=[],nGroup=10):
	'''relative conductivities'''
	#nGroup=10
	#Names = genfromtxt('grads_hdr.csv', delimiter=' ',dtype='string').tolist()
	#volInt = genfromtxt('grads.csv', delimiter=',',dtype='double',skip_header=1)[1:-1]
	#SinglePhaseRes=[[]]
	#SinglePhaseRes[0] = genfromtxt('grads_SP.csv', delimiter=',',dtype='double')

	# compute integrals
	t=groupTime(volInt,0,nGroup)
	alpha1=groupTimeAvg(volInt,Names,('alpha'),inds,nGroup)
	Pc=groupTimeAvg(volInt,Names,('Pc'),inds,nGroup)
	UzAvg=groupTimeAvg(volInt,Names,('U'),inds,nGroup)
	f_1=groupTimeAvg(volInt,Names,('f_1'),inds,nGroup)
	Uz1Avg=(UzAvg*f_1)
	Uz2Avg=UzAvg - Uz1Avg
	viscEV=groupTimeAvg(volInt,Names,('viscE'),inds,nGroup)
	viscEV1=groupTimeAvg(volInt,Names,('viscE_1'),inds,nGroup)
	viscEV2=viscEV - viscEV1 - 1e-6
	UzAvgSP=groupTimeAvg(SinglePhaseRes,Names,'U',inds,1)
	viscEVSP=groupTimeAvg(SinglePhaseRes,Names,(('viscE')),inds,1)

	#viscosity
	mu1=viscEV2*0.0+1.0
	mu2=mu1
	muSP=mu1
	try : mu1=groupTimeAvg(volInt,Names,(('mu1')),inds,nGroup)
	except : pass
	try : mu2=groupTimeAvg(volInt,Names,(('mu2')),inds,nGroup)
	except : pass
	try : muSP=groupTimeAvg(SinglePhaseRes,Names,(('mu1')),inds,1)
	except : pass
	print(("visc:",mu1[-1],mu2[-1],muSP[-1]))

	# compute and return Sw, Pc, Krw and Kro
	viscEVSP=viscEVSP[-1]/muSP[-1]
	SwPcKrs=[[],[],[],[],[]]
	SwPcKrs[0]=1.-alpha1
	SwPcKrs[1]=Pc
	SwPcKrs[3]=(Uz1Avg*Uz1Avg) / (viscEV1/mu1-1e-16) / ((UzAvgSP*UzAvgSP) / viscEVSP)
	SwPcKrs[2]=(Uz2Avg*Uz2Avg) / (viscEV2/mu2-1e-16) / ((UzAvgSP*UzAvgSP) / viscEVSP)
	SwPcKrs[4]=t
	return SwPcKrs




def plotMisc(xxx=None,indx=None,delNGroup=None):
	Colours=['-.1x', ':4+', '--3o', '-5*', '-2',   '--', '2+', '3x', '4+', '5*',   '1o', '2x', '3+', '4*', '5o',   '1x', '2+', '3*', '4o', '5x', ]

	#################################################   111111111 TOTAL FORCE ########################################################
	if False:
		viscFzV=groupTimeAvg(volInt,Names,('viscz'),inds,nGroup)
		viscFzV1=groupTimeAvg(volInt,Names,('viscz_1'),inds,nGroup)
		viscFzV2=viscFzV - viscFzV1

		viscInterf_1=groupTimeAvg(volInt,Names,('viscInterf_1'),inds,nGroup)

		dPdzV=groupTimeAvg(volInt,Names,('dpddz'),inds,nGroup)
		dPdzV1=groupTimeAvg(volInt,Names,('dpddz_1'),inds,nGroup)
		dPdzV2=dPdzV - dPdzV1

		dPczV=groupTimeAvg(volInt,Names,('dpcdz'),inds,nGroup)
		dPczV1=groupTimeAvg(volInt,Names,('dpcdz_1'),inds,nGroup)
		dPczV2=dPczV - dPczV1

		# phiUxS =-1000*surfInt(1:end,SrfI('phiU:0')) ;
		# phiUyS =-1000*surfInt(1:end,SrfI('phiU:1')) ;
		phiUzS=groupTimeAvg(volInt,Names,('phiu'),inds,nGroup)

		phiUzS1=groupTimeAvg(volInt,Names,('phiu_1'),inds,nGroup)

		phiUzS1=phiUzS - phiUzS1

		# dlmwrite(["plotTotalForce"  ".csv"], t', '-append') ;
		# dlmwrite(["plotTotalForce"  ".csv"], viscFzV', '-append') ;
		# dlmwrite(["plotTotalForce"  ".csv"], viscFzV'./UzAvg', '-append') ;
		# dlmwrite(["plotTotalForce"  ".csv"], dPdzV'./UzAvg', '-append') ;
		hold('off')
		plot(dot(t,10**3),(multiply(dpEd,Volume) / ((qBack))) / Pc_ref,Colours[2,:],'LineWidth',2,'markersize',2)
		hold('on')
		plot(dot(t,10**3),(multiply(dpEc,Volume) / ((qBack))) / Pc_ref,Colours[1,:],'LineWidth',2,'markersize',2)
		plot(dot(t,10**3),(multiply(viscEV,Volume) / ((qBack))) / Pc_ref,Colours[3,:],'LineWidth',2,'markersize',2)
		# plot(t,(dPdzV2+dPczV2)/10000000.0,Colours(3,:),'LineWidth',2);
		# plot(t,(dPczV1+dPczV2)/10000000.0,Colours(2,:),'LineWidth',2);
		# plot(t*10^3,((viscEV+dpEc+dpEd+phiE).*Volume./((UzAvg.*AreaCross)))./Pc_ref,Colours(6,:),'LineWidth',2,'markersize',2);
		# max(viscFzV)
		plot(dot(t,10**3),1 - dot(alpha1,1.0),Colours[5,:],'LineWidth',2)
		disp('________________TOTAL______________-')
		# ylim([-2,2])
		# yMin=0;yMax=-100000;
		# meanTmp=mean(dPdzV1/10000.0);	 yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		# meanTmp=mean(dPczV1/10000.0);	 yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		# meanTmp=mean(dPdzV2/10000.0);	 yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		# meanTmp=mean(dPczV2/10000.0);	 yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		#ufactor=yMax/mean(UzAvg./( dPdzV+dPczV)*dppq_sp(iCV));
		# yMax=ceil(yMax*10)/10;
		# yMin=floor(10*yMin)/10;
		#afactor=yMax;
		# delPcelZ =-groupTimeAvg(volInt,Names,(["-delPcelZ"]),inds, nGroup  ) ;
		# plot(t, (viscFzV1+dPdzV1+dPczV1)/10000.0,Colours(5,:),'LineWidth',2);

		# plot(t, (viscFzV2+dPdzV2+dPczV2)/10000.0,Colours(6,:),'LineWidth',2);
		# xlim([0 10000])
		# ylim([-0.03 0.06])
		# set (gca, "interpreter", "latex")
		# set(gca,'FontSize',14);
		# set(gca,'FontName',"LMRoman10");
		set(gca,'xminortick','on','yminortick','on','xgrid','on','ygrid','on')
		# xlim([0.00,0.01])
		# yLims=ylim();  ylim(yLims*0.5)
		# yScaleSUGGEST=[3*yMin 3*yMax]
		# ylim(yScale(3,:))
		xlabel('t  ({m}s)')
		# legend(fileNames,"Location","North");
		Leg=legend('{ P_{d}} (kPa)','{ P_{c}} (kPa)','{ P_{ mu}} (kPa)','{ dE_{k;Adv}|Qdt} (kPa)','{ S_{w}}  ','Location','southeast')

		# set(FL1,'Spacing',2);
		# set (h, "interpreter", "tex")
		# legend boxon
		# set(gca,"keypos", 2)
		# # title (str(iCV));
		# print(["plotTotalForce" str(iCV) ".png"],"-dpng","-r90","-FLMRoman10:14")  ;
		# print(["plotTotalForce" str(iCV) ".eps"], "-deps", "-color","-FLMRoman10:14")  ;
		#prin t2pdflatex('1','plotTotalForce',str(iCV))
		hold('off')
	#################################################   222222222  EACH PHASE ENERGY LOSS 22222222 ########################################################
	if False:
		plot(dot(t,10 ** 3),(multiply(multiply(- (f_1 / (abs(f_1) + 1e-06)),(viscEV1)),Volume) / (multiply((Uz1Avg + 1e-12),AreaCross))) / Pc_ref,Colours[1,:],'LineWidth',2,'markersize',2)
		hold('on')
		plot(dot(t,10 ** 3),(multiply(multiply(- (f_2 / (abs(f_2) + 1e-06)),(viscEV2)),Volume) / (multiply((Uz2Avg + 1e-12),AreaCross))) / Pc_ref,Colours[3,:],'LineWidth',2,'markersize',2)
		# plot(t*10^3,( dpEc2.*Volume./(Uz2Avg.*AreaCross))./Pc_ref,Colours(4,:),'LineWidth',2,'markersize',2);
		# max(viscFzV)
		# plot(t,dPczV2/10000000.0,Colours(4,:),'LineWidth',2);
		plot(cat(0),cat(-1),Colours[5,:],'LineWidth',4)
		#plot(t*10^3,( -(f_2./(abs(f_2)+1e-6)).*(dpEd2+dpEc2).*Volume./((Uz2Avg+1e-12).*AreaCross))./Pc_ref,Colours(3,:),'LineWidth',1,'markersize',0);
		disp('________________EACH______________-')
		# ylim([0 10])
		# xlim([0.00,900])
		# plot(t*10^3, alpha1*0,'-0','LineWidth',1);
		set(gca,'xminortick','on','yminortick','on','xgrid','on','ygrid','on')
		# yLims=ylim()
		# ylim(yLims*0.5)
		# yScaleSUGGEST=[3*yMin 3*yMax]
		# ylim(yScale(3,:))
		xlabel('{ t} (ms)')
		ylabel('{  P_{}} (kPa)')
		Leg=legend('{ P_{,nw}} ','{ P_{,w}} ','{ S_{w}}  ','Location','northeast')

		haxes1=copy(gca)

		haxes1_pos=get(gca,'Position')

		haxes2=axes('Position',haxes1_pos,'xaxisLocation','top','yaxisLocation','right','Color','none')

		plot(dot(t,10 ** 3),dot((1 - alpha1),1),'Parent',haxes2,Colours[5,:],'LineWidth',4)
		#set(haxes2,'xlim',[1 300]);
		# set(haxes2,'xlabel'," ");
		#					 xmin ymin xmax ymax
		set(haxes1,'position',cat(7,3,3,- 2))
		set(haxes2,'position',cat(7,3,3,- 2))
		set(haxes2,'yaxisLocation','right')
		set(haxes2,'XTickLabel','')
		# set (haxes1,'defaultaxesposition', [0.05, 0.1, 0.9, 0.85])
		# set (haxes2,'defaultaxesposition', [0.05, 0.1, 0.9, 0.85])
		set(haxes2,'ylabel','{ S_{w}}  ')
		# print(["plotEachForce" str(iCV) ".png"],"-dpng","-r90","-FLMRoman10:14");
		# print(["plotEachForce" str(iCV) ".eps"], "-deps", "-color","-FLMRoman10:14")  ;
		#prin t2pdflatex('2',cat('plotEachForce',str(iCV)))
		hold('off')
		# continue
		clf
	#################################################   222222222  EACH PHASEnoxs1 22222222 ########################################################
	if False:
		plot(dot(t,10 ** 3),(multiply(dpEd1,Volume) / (multiply(Uz1Avg,AreaCross))) / Pc_ref,Colours[1,:],'LineWidth',2,'markersize',2)
		hold('on')
		plot(dot(t,10 ** 3),(multiply(dpEd2,Volume) / (multiply(Uz2Avg,AreaCross))) / Pc_ref,Colours[2,:],'LineWidth',2,'markersize',2)
		plot(dot(t,10 ** 3),(multiply(dpEc1,Volume) / (multiply(Uz1Avg,AreaCross))) / Pc_ref,Colours[3,:],'LineWidth',2,'markersize',2)
		plot(dot(t,10 ** 3),(multiply(dpEc2,Volume) / (multiply(Uz2Avg,AreaCross))) / Pc_ref,Colours[4,:],'LineWidth',2,'markersize',2)
		# plot(t,dPczV2/10000000.0,Colours(4,:),'LineWidth',2);
		plot(dot(t,10 ** 3),1 - dot(alpha1,1),Colours[5,:],'LineWidth',4)
		disp('________________EACH______________-')
		# ylim([-0.5 2])
		# xlim([0.00,1000])
		plot(dot(t,10 ** 3),dot(alpha1,0),'-0','LineWidth',1)
		set(gca,'xminortick','on','yminortick','on','xgrid','on','ygrid','on')
		# yLims=ylim()
		# ylim(yLims*0.5)
		# yScaleSUGGEST=[3*yMin 3*yMax]
		# ylim(yScale(3,:))
		xlabel('{ t} ({m}s)')
		# legend(fileNames,"Location","North");
		Leg=legend('{  P_{d,nw}} (kPa)','{  P_{d,w}} (kPa)','{  P_{c,nw}} (kPa)','{  P_{c,w}} (kPa)','{ S_{w}}  ','Location','northeast')

		# pr int(["plotEachForce" str(iCV) ".eps"], "-deps", "-color","-FLMRoman10:14")  ;
		#prin t2pdflatex('2',('plotEachForceNoXf'+str(iCV)))
		hold('off')

	#################################################   333333333  REL PermNoPc ########################################################
	if False:
		UzAvgSP=groupTimeAvg(SinglePhaseRes,Names,(cat('U')),inds,1)
		viscEVSP=groupTimeAvg(SinglePhaseRes,Names,(cat('viscE')),inds,1)
		UzAvgSP=copy(UzAvgSP)
		viscEVSP=viscEVSP[end()]
		plot(1.0 - alpha1,multiply(Uz1Avg,UzAvg) / (viscEV) / (dot(UzAvgSP,UzAvgSP) / viscEVSP),Colours[1,:],'LineWidth',2,'markersize',3)
		hold('on')
		plot(1.0 - alpha1,multiply(Uz2Avg,UzAvg) / (viscEV) / (dot(UzAvgSP,UzAvgSP) / viscEVSP),Colours[3,:],'LineWidth',2,'markersize',3)
		plot(1.0 - alpha1,f_1,Colours[2,:],'LineWidth',2,'markersize',1)
		disp('__________pppppppppppp______________-')
		set(gca,'xminortick','on','yminortick','on','xgrid','on','ygrid','on')
		xlim(cat(- 0.0001,1.0001))
		ylim(cat(0.0001,1.1))
		xlabel('{ S_w}')
		legend('  {   k_{r,nw}}','  {   k_{r,w}}','  {   f_{nw}}','Location','west')
		hold('off')
	#################################################   222 VISC  REL PERM ########################################################






if __name__ == "__main__":
	caseName='.'
	nCntrlVols=1
	if os.path.exists(caseName+'/grads.csv') :
	disp(caseName)
	Names = genfromtxt(caseName+'/grads_hdr.csv', delimiter=' ',dtype='string').tolist()
	volInt = genfromtxt(caseName+'/grads.csv', delimiter=',',dtype='double',skip_header=1)[1:-1]
	SinglePhaseRes=[[]]
	SinglePhaseRes[0] = genfromtxt(caseName+'/grads_SP.csv', delimiter=',',dtype='double')[-1]
	relPrems = gradsToRelPermSwPcKrwKroTRev(Names,volInt,SinglePhaseRes,np.arange(1,nCntrlVols),10) #cas1Ncors[i][k]/2+
	relPrems[1][:] *= 1/0.03
	np.savetxt(caseName+'_relPerms.tsv', relPrems, delimiter='\t')

