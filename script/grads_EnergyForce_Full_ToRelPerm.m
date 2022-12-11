% outdated! usegroupGrads.py instead

function Error=processPlotRelPerm();
% yScale
format compact
     hold off


% dx=0.000012/20;
Colours=[ ...
'-.1x'; ':4+'; '--3o'; '-5*'; '-2'; ...
'--'; '2+'; '3x'; '4+'; '5*'; ...
'1o'; '2x'; '3+'; '4*'; '5o'; ...
'1x'; '2+'; '3*'; '4o'; '5x'; ];
% Colours= cellstr (Colours_);
% Colours(1)
% Colours_(1,:)
% fileNames{1}
% length(fileNames)

Pc_ref=1000;

nGroup=10



volInt=csvread('grads.csv');
Names = cellstr(textread ('grads_hdr.csv','%s'));
SinglePhaseRes=csvread('grads_SP.csv');
SinglePhaseRes=SinglePhaseRes(end,:);


% i=floor((length(Names)-10)/30)+1 ;     i=i-1
% while(i>  5   )
% i=i-1
ihas=[1 2 3 4 5 6 7 8];
for iiCV=1:length(ihas)
	iha=ihas(iiCV);
	iCV=ihas(iiCV);

	t=groupTime(volInt,[1], nGroup  ) ;
	Volume=groupTimeSum(volInt,Names,(['vol']),iha, nGroup  ) ;
	alpha1=groupTimeAvg(volInt,Names,(['alpha']),iha, nGroup  ) ;
	alpha2=1.-alpha1;
	Volume1=Volume.*alpha1;
	Volume2=Volume-Volume1;

	z1=groupTimeMin(volInt,Names,(['x1']),iha, nGroup  ) ;
	z2=groupTimeMax(volInt,Names,(['x2']),iha, nGroup  ) ;
	xDropAvg=groupTimeAvg(volInt,Names,(['xDropAvg']),iha, nGroup  ) ;
	xDrop1=groupTimeMin(volInt,Names,(['xDrop1']),iha, nGroup  ) ;
	xDrop2=groupTimeMax(volInt,Names,(['xDrop2']),iha, nGroup  ) ;


	Volume=Volume+1e-18;
	Volume1=Volume1+1e-18;
	Volume2=Volume2+1e-18;



	% AreaCross=groupTime(volInt,VolI(Names,'Volume'))/DeltaX;

	% ViscCorrection=1.+0.*(1.-0.5*(dPdzV1+dPczV1)./groupTime(volInt,VolI(Names,'Volume'))*dx/2 ...
	% .*((surfInt(1:end,SrfI('Area'))-2*AreaCross)*0.66)./viscFzV);



	% UzAvg=groupTime(volInt,VolI(Names,'U:2'));
	% UzTubeV1=-(AreaCross/(2*pi()).*( viscFzV.*ViscCorrection)/0.001);% + dPdzV1+dPczV1
	% UzTubeV=-0.4*(AreaCross./(2*pi()).*( -dPdzV-dPczV)/0.001);% + dPdzV1+dPczV1
	% UzTubeV1=-(AreaCross./(2*pi()).*( -dPdzV1-dPczV1).*Volume./Volume1/0.001);% + dPdzV1+dPczV1
	% UzTubeV2=-(AreaCross./(2*pi()).*( -dPdzV2-dPczV2).*Volume./Volume2/0.001);% + dPdzV1+dPczV1


	UzAvg=groupTimeAvg(volInt,Names,(['U']),iha, nGroup  ) ;
	f_1=groupTimeAvg(volInt,Names,(['f_1']),iha, nGroup  ) ;
	f_2=1.-f_1;
	Uz1Avg=UzAvg.*f_1;
	Uz2Avg=UzAvg-Uz1Avg;
	% QPQTube=UzAvg./( dPdzV+dPczV)*dppq_sp(iCV);



	% z2-z1
	AreaCross=Volume./(z2-z1);
	AreaCross1=Volume1./(z2-z1);
	AreaCross2=Volume2./(z2-z1);




	qBack=-groupTime(volInt,VolI(Names,['QIn']), nGroup  ) ;
	qFront=groupTime(volInt,VolI(Names,['QOut']), nGroup  ) ;

	porevolume=[0.5*(t(2:end)+t(1:end-1));cumsum(0.5*(qBack(2:end)+qBack(1:end-1)).*diff(t)./Volume(2:end));alpha1(2:end)];



	viscEV=groupTimeAvg(volInt,Names,(['viscE']),iha, nGroup  ) ;
	viscEV1=groupTimeAvg(volInt,Names,(['viscE_1']),iha, nGroup  )
	viscEV2=viscEV-viscEV1-1e-12
	dpEc =groupTimeAvg(volInt,Names,(['dpEc']),iha, nGroup  ) ;
	dpEc1 =groupTimeAvg(volInt,Names,(['dpEc_1']),iha, nGroup  ) ;
	dpEc2 =dpEc-dpEc1;
	dpEd =groupTimeAvg(volInt,Names,(['dpEd']),iha, nGroup  ) ;
	dpEd1 =groupTimeAvg(volInt,Names,(['dpEd_1']),iha, nGroup  ) ;
	dpEd2 =dpEd-dpEd1;
	phiE =groupTimeAvg(volInt,Names,(['phiE']),iha, nGroup  ) ;
	phiE1 =groupTimeAvg(volInt,Names,(['phiE_1']),iha, nGroup  ) ;
	phiE2 =dpEd-dpEd1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   111111111 TOTAL FORCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if (0)

		viscFzV=groupTimeAvg(volInt,Names,(['viscz']),iha, nGroup  ) ;
		viscFzV1=groupTimeAvg(volInt,Names,(['viscz_1']),iha, nGroup  ) ;
		viscFzV2=viscFzV-viscFzV1 ;
		viscInterf_1=groupTimeAvg(volInt,Names,(['viscInterf_1']),iha, nGroup  ) ;

		% dPdxV =-groupTime(volInt,VolI(Names,'gPd:0'));
		% dPdyV =-groupTime(volInt,VolI(Names,'gPd:1'));
		dPdzV =groupTimeAvg(volInt,Names,(['dpddz']),iha, nGroup  ) ;
		dPdzV1 =groupTimeAvg(volInt,Names,(['dpddz_1']),iha, nGroup  ) ;
		dPdzV2 =dPdzV-dPdzV1;
		% dPdzV1 =-surfInt(1:end,SrfI('pxN:2'));

		% dPcxV =groupTime(volInt,VolI(Names,'gPc:0'));
		% dPcyV =groupTime(volInt,VolI(Names,'gPc:1'));
		dPczV =groupTimeAvg(volInt,Names,(['dpcdz']),iha, nGroup  ) ;
		dPczV1 =groupTimeAvg(volInt,Names,(['dpcdz_1']),iha, nGroup  ) ;
		dPczV2 =dPczV-dPczV1;
		% dPczV1 =-surfInt(1:end,SrfI('PcxN:2'));

		% phiUxS =-1000*surfInt(1:end,SrfI('phiU:0')) ;
		% phiUyS =-1000*surfInt(1:end,SrfI('phiU:1')) ;
		phiUzS =groupTimeAvg(volInt,Names,(['phiu']),iha, nGroup  ) ;
		phiUzS1 =groupTimeAvg(volInt,Names,(['phiu_1']),iha, nGroup  ) ;
		phiUzS1 =phiUzS-phiUzS1;


		% dlmwrite(['plotTotalForce" ".csv"], ["CV" num2str(iCV) ":"], '-append') ;
		% dlmwrite(["plotTotalForce"  ".csv"], t', '-append') ;
		% dlmwrite(["plotTotalForce"  ".csv"], viscFzV', '-append') ;
		% dlmwrite(["plotTotalForce"  ".csv"], viscFzV'./UzAvg', '-append') ;
		% dlmwrite(["plotTotalForce"  ".csv"], dPdzV'./UzAvg', '-append') ;
		hold off
		plot(t*10^3,(dpEd.*Volume./((qBack)))./Pc_ref,Colours(2,:),'LineWidth',2,'markersize',2);
		hold on
		plot(t*10^3,(dpEc.*Volume./((qBack)))./Pc_ref,Colours(1,:),'LineWidth',2,'markersize',2);
		plot(t*10^3,(viscEV.*Volume./((qBack)))./Pc_ref,Colours(3,:),'LineWidth',2,'markersize',2);
		% plot(t*10^3,(phiE.*Volume./((qBack)))./Pc_ref,Colours(4,:),'LineWidth',2,'markersize',2);
		% plot(t,(dPdzV2+dPczV2)/10000000.,Colours(3,:),'LineWidth',2);
		% plot(t,(dPczV1+dPczV2)/10000000.,Colours(2,:),'LineWidth',2);
		% plot(t*10^3,((viscEV+dpEc+dpEd+phiE).*Volume./((UzAvg.*AreaCross)))./Pc_ref,Colours(6,:),'LineWidth',2,'markersize',2);
		% max(viscFzV)
		plot(t*10^3, 1-alpha1*1.,Colours(5,:),'LineWidth',2);

		disp('________________TOTAL______________-');


		% xlim([0.00,900])
		% ylim([-2,2])

		% yMin=0;yMax=-100000;
		% meanTmp=mean(dPdzV1/10000.);     yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		% meanTmp=mean(dPczV1/10000.);     yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		% meanTmp=mean(dPdzV2/10000.);     yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		% meanTmp=mean(dPczV2/10000.);     yMin=min(meanTmp,yMin); yMax=max(yMax,meanTmp);
		%ufactor=yMax/mean(UzAvg./( dPdzV+dPczV)*dppq_sp(iCV));
		% yMax=ceil(yMax*10)/10;
		% yMin=floor(10*yMin)/10;
		%afactor=yMax;

		% delPcelZ =-groupTimeAvg(volInt,Names,(["-delPcelZ"]),iha, nGroup  ) ;
		% plot(t, (viscFzV1+dPdzV1+dPczV1)/10000.,Colours(5,:),'LineWidth',2);
		%
		% plot(t, (viscFzV2+dPdzV2+dPczV2)/10000.,Colours(6,:),'LineWidth',2);

		% xlim([0 10000])
		% ylim([-0.03 0.06])

		% set (gca, "interpreter", "latex")
		% set(gca,'FontSize',14);
		% set(gca,'FontName',"LMRoman10");

		set (gca, "xminortick", "on", "yminortick", "on", "xgrid", "on", "ygrid", "on") ;
		% axis([0 8e-11 -4000  4000 ]);
		% xlim([0.00,0.01])
		% yLims=ylim()
		% ylim(yLims*0.5)
		% yScaleSUGGEST=[3*yMin 3*yMax]
		% ylim(yScale(3,:))
		xlabel("t  ({m}s)");
		% ylabel("Force (10^6N/m^3)");
		% legend(fileNames,"Location","North");
		Leg=legend(...
		'{/LMRomanUnsl10 𝛥P_{d}} (kPa)',...
		'{/LMRomanUnsl10 𝛥P_{c}} (kPa)',...
		'{/LMRomanUnsl10 𝛥P_{\mu}} (kPa)',...
		'{/LMRomanUnsl10 dE_{k;Adv}|Qdt} (kPa)',...
		% % '{/LMRomanUnsl10 dE_{k,Acc.}/Qdt} (kPa)',...
		'{/LMRomanUnsl10 S_{w}}  ',...
		% % 'Visc 2\n',...
		'Location','southeast');


		% FL1= findall(Leg,'-property','Spacing');
		% set(FL1,'Spacing',2);
		% set (h, "interpreter", "tex")
		% legend boxon
		% set(gca,"keypos", 2)
		% % title (num2str(iCV));

		% print(["plotTotalForce" num2str(iCV) ".png"],"-dpng","-r90","-FLMRoman10:14")  ;
		% print(["plotTotalForce" num2str(iCV) ".eps"], "-deps", "-color","-FLMRoman10:14")  ;
		print2pdflatex("1",["plotTotalForce" num2str(iCV)]);

		hold off

	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   222222222  EACH PHASE ENERGY LOSS 22222222 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if (0)
		plot(t*10^3,( -(f_1./(abs(f_1)+1e-6)).*(viscEV1).*Volume./((Uz1Avg+1e-12).*AreaCross))./Pc_ref,Colours(1,:),'LineWidth',2,'markersize',2);
		hold on

		plot(t*10^3,( -(f_2./(abs(f_2)+1e-6)).*(viscEV2).*Volume./((Uz2Avg+1e-12).*AreaCross))./Pc_ref,Colours(3,:),'LineWidth',2,'markersize',2);
		% plot(t*10^3,( dpEc1.*Volume./(Uz1Avg.*AreaCross))./Pc_ref,Colours(3,:),'LineWidth',2,'markersize',2);
		% plot(t*10^3,( dpEc2.*Volume./(Uz2Avg.*AreaCross))./Pc_ref,Colours(4,:),'LineWidth',2,'markersize',2);
		% max(viscFzV)
		% plot(t,dPczV2/10000000.,Colours(4,:),'LineWidth',2);

		plot([0], [-1], Colours(5,:),'LineWidth',4);

		%plot(t*10^3,( -(f_1./(abs(f_1)+1e-6)).*(dpEd1+dpEc1).*Volume./((Uz1Avg+1e-12).*AreaCross))./Pc_ref,Colours(1,:),'LineWidth',1,'markersize',0);
		%plot(t*10^3,( -(f_2./(abs(f_2)+1e-6)).*(dpEd2+dpEc2).*Volume./((Uz2Avg+1e-12).*AreaCross))./Pc_ref,Colours(3,:),'LineWidth',1,'markersize',0);

		disp("________________EACH______________-");

		% xlim([1 80])
		% ylim([0 10])
		% xlim([0.00,900])
		% plot(t*10^3, alpha1*0,'-0','LineWidth',1);

		set (gca, "xminortick", "on", "yminortick", "on", "xgrid", "on", "ygrid", "on") ;
		% axis([0 8e-11 -4000  4000 ]);
		% yLims=ylim()
		% ylim(yLims*0.5)
		% yScaleSUGGEST=[3*yMin 3*yMax]
		% ylim(yScale(3,:))
		xlabel("{/LMRomanUnsl10 t} (ms)");
		ylabel("{/LMRomanUnsl10  𝛥P_{μ}} (kPa)");
		% legend(fileNames,"Location","North");

		Leg=legend(...
		"{/LMRomanUnsl10 𝛥P_{μ,nw}} ",... %𝜇
		"{/LMRomanUnsl10 𝛥P_{μ,w}} ",... %𝜇
		% "\ \ {/LMRomanUnsl10  f_{nw}𝛥P_{c,nw}} (kPa)",...
		% "\ \ {/LMRomanUnsl10  f_{w}𝛥P_{c,w}} (kPa)",...
				"{/LMRomanUnsl10 S_{w}}  ",...
		"Location","northeast");


		haxes1 = gca; % handle to axes
		haxes1_pos = get(gca,'Position'); % store position of first axes

		haxes2 = axes('Position',haxes1_pos,...
		  'xaxisLocation','top',...
			  'yaxisLocation','right',...
			  'Color','none');
		plot(t*10^3, (1-alpha1)*1,'Parent',haxes2, Colours(5,:),'LineWidth',4);

		%set(haxes2,'ylim',[. 1]);
		%set(haxes2,'xlim',[1 300]);
		% set(haxes2,'xlabel'," ");
		%                     xmin ymin xmax ymax
		set (haxes1, "position",  [7, 3, 3, -2]);
		set (haxes2, "position",  [7, 3, 3, -2]);
		set(haxes2,'yaxisLocation','right')
		set(haxes2,'XTickLabel','')
		% set (haxes1,'defaultaxesposition', [0.05, 0.1, 0.9, 0.85])
		% set (haxes2,'defaultaxesposition', [0.05, 0.1, 0.9, 0.85])
		set (haxes2,"ylabel", "{/LMRomanUnsl10 S_{w}}  ")



		% print(["plotEachForce" num2str(iCV) ".png"],"-dpng","-r90","-FLMRoman10:14");
		% print(["plotEachForce" num2str(iCV) ".eps"], "-deps", "-color","-FLMRoman10:14")  ;
		print2pdflatex("2",["plotEachForce" num2str(iCV)]);

		hold off

		% continue
		clf


	end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   222222222  EACH PHASEnoxs1 22222222 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (0)

		plot(t*10^3,( dpEd1.*Volume./(Uz1Avg.*AreaCross))./Pc_ref,Colours(1,:),'LineWidth',2,'markersize',2);
		hold on
		plot(t*10^3,( dpEd2.*Volume./(Uz2Avg.*AreaCross))./Pc_ref,Colours(2,:),'LineWidth',2,'markersize',2);
		plot(t*10^3,( dpEc1.*Volume./(Uz1Avg.*AreaCross))./Pc_ref,Colours(3,:),'LineWidth',2,'markersize',2);
		plot(t*10^3,( dpEc2.*Volume./(Uz2Avg.*AreaCross))./Pc_ref,Colours(4,:),'LineWidth',2,'markersize',2);
		% max(viscFzV)
		% plot(t,dPczV2/10000000.,Colours(4,:),'LineWidth',2);
		plot(t*10^3, 1-alpha1*1,Colours(5,:),'LineWidth',4);

		disp("________________EACH______________-");

		% xlim([0 10000])
		% ylim([-0.5 2])
		 % xlim([0.00,1000])
		  plot(t*10^3, alpha1*0,'-0','LineWidth',1);


		  set (gca, "xminortick", "on", "yminortick", "on", "xgrid", "on", "ygrid", "on") ;
		  % axis([0 8e-11 -4000  4000 ]);
		 % yLims=ylim()
		 % ylim(yLims*0.5)
		 % yScaleSUGGEST=[3*yMin 3*yMax]
		 % ylim(yScale(3,:))
		  xlabel("{/LMRomanUnsl10 t} ({m}s)");
		  % ylabel("Force (10^4N/m^3)");
		  % legend(fileNames,"Location","North");
		  Leg=legend(...
		  "\ {/LMRomanUnsl10  𝛥P_{d,nw}} (kPa)",... %𝜇
		  "{/LMRomanUnsl10  𝛥P_{d,w}} (kPa)",... %𝜇
			 "{/LMRomanUnsl10  𝛥P_{c,nw}} (kPa)",...
			 "{/LMRomanUnsl10  𝛥P_{c,w}} (kPa)",...
						"{/LMRomanUnsl10 S_{w}}  ",...
		  "Location","northeast");


			% print(["plotEachForce" num2str(iCV) ".png"],"-dpng","-r90","-FLMRoman10:14");
			% print(["plotEachForce" num2str(iCV) ".eps"], "-deps", "-color","-FLMRoman10:14")  ;
		  print2pdflatex("2",["plotEachForceNoXf" num2str(iCV)]);

		  hold off

	% continue

	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   333333333  REL PERM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
	% UzAvgSP=SinglePhaseRes(end,VolI(Names,["S" num2str(iCV) "-U"]))
	% viscEVSP=SinglePhaseRes(end,VolI(Names,["S" num2str(iCV) "-viscE"]))
	UzAvgSP=groupTimeAvg(SinglePhaseRes,Names,(['U']),iha,1);
	viscEVSP=groupTimeAvg(SinglePhaseRes,Names,(['viscE']),iha,1);
	dpEdVSP=groupTimeAvg(SinglePhaseRes,Names,(['dpEd']),iha,1);


	viscEVSP=viscEVSP(length(viscEVSP))



	"RelPerm"
	[1.-alpha1;   Uz1Avg.*Uz1Avg./( viscEV1) /(UzAvgSP*UzAvgSP/viscEVSP);   Uz2Avg.*Uz2Avg./( viscEV2) /(UzAvgSP*UzAvgSP/viscEVSP)]'

	semilogy(1.-alpha1,Uz1Avg.*Uz1Avg./( viscEV1) /(UzAvgSP*UzAvgSP/viscEVSP),Colours(1,:),'LineWidth',2,'markersize',3);
	hold on
	semilogy(1.-alpha1,Uz2Avg.*Uz2Avg./( viscEV2)   /(UzAvgSP*UzAvgSP/viscEVSP),Colours(3,:),'LineWidth',2,'markersize',3);
	semilogy(1.-alpha1,t,Colours(2,:),'LineWidth',2,'markersize',0);

	disp('__________pppppppppppp______________-');


	dlmwrite(['relPerm', num2str(iCV), '.csv'], [1.-alpha1;Uz1Avg.*Uz1Avg./( viscEV1)   /(UzAvgSP*UzAvgSP/viscEVSP) ...
	;Uz2Avg.*Uz2Avg./( viscEV2)   /(UzAvgSP*UzAvgSP/viscEVSP)]');


	%plot(1.-alpha1,Uz1Avg.*Uz1Avg./( dpEc1+dpEd1)   /(UzAvgSP*UzAvgSP/dpEdVSP),Colours(1,:),'LineWidth',1,'markersize',0);
	%plot(1.-alpha1,Uz2Avg.*Uz2Avg./( dpEc2+dpEd2)   /(UzAvgSP*UzAvgSP/dpEdVSP),Colours(3,:),'LineWidth',1,'markersize',0);

	set (gca, 'xminortick', 'on', 'yminortick', 'on', "xgrid", "on", "ygrid", "on") ;
	% axis([0 0.02 -0.2  1.2 ]);
	xlim([-0.0001,1.0001]);
	% yLims=ylim()
	% ylim(yLims*0.5)
	ylim([0.001,1.000]);
	xlabel("{/LMRomanUnsl10 S_w}");
	% ylabel("\{ \/LMRomanUnsl10 q/𝛥P / (q/𝛥P)_{SP} \}");
	% legend(fileNames,'Location','North');
	legend(...
	"  k_{r,o}",...
	"  k_{r,w}",...
	"  t (s)",...
	"Location","west");

	print2pdflatex("3",["plotRelPerm" num2str(iCV)]);

	hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   333333333  REL PermNoPc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if (0)
		% UzAvgSP=SinglePhaseRes(end,VolI(Names,["S" num2str(iCV) "-U"]))
		% viscEVSP=SinglePhaseRes(end,VolI(Names,["S" num2str(iCV) "-viscE"]))
		UzAvgSP=groupTimeAvg(SinglePhaseRes,Names,(["U"]),iha,1);
		viscEVSP=groupTimeAvg(SinglePhaseRes,Names,(["viscE"]),iha,1);
		UzAvgSP=UzAvgSP
		viscEVSP=viscEVSP(end)

		plot(1.-alpha1,Uz1Avg.*UzAvg./( viscEV)   /(UzAvgSP*UzAvgSP/viscEVSP),Colours(1,:),'LineWidth',2,'markersize',3);
		hold on
		plot(1.-alpha1,Uz2Avg.*UzAvg./( viscEV)   /(UzAvgSP*UzAvgSP/viscEVSP),Colours(3,:),'LineWidth',2,'markersize',3);
		plot(1.-alpha1,f_1,Colours(2,:),'LineWidth',2,'markersize',1);
		% plot(alpha1,10*Uz2Avg./( viscEV2 )    /(UzAvgSP/viscEVSP),Colours(4,:),'LineWidth',2,'markersize',2);

		disp('__________pppppppppppp______________-');

		% dPdzV1+dPczV1+viscFzV.*ViscCorrection+

		set (gca, 'xminortick', 'on', 'yminortick', 'on', "xgrid", "on", "ygrid", "on") ;
		% axis([0 0.02 -0.2  1.2 ]);
		xlim([-0.0001,1.0001]);
		% yLims=ylim()
		% ylim(yLims*0.5)
		ylim([0.0001,1.1]);
		xlabel("{/LMRomanUnsl10 S_w}");
		% ylabel("\{ \/LMRomanUnsl10 q/𝛥P / (q/𝛥P)_{SP} \}");
		% legend(fileNames,'Location','North');
		legend(...
		"  {/LMRomanUnsl10   k_{r,nw}}",...
		"  {/LMRomanUnsl10   k_{r,w}}",...
		"  {/LMRomanUnsl10   f_{nw}}",...
		"Location","west");
		% title (num2str(iCV));

		print2pdflatex("3",["plotRelPermNoPc" num2str(iCV)]);

		hold off
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   222 VISC  REL PERM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














	end
end







% function  returnedDataInd=    SrfI(Name)
% Names={"Time", ...
% "gUyxN1","phiU1:0","phiU1:1","phiU1:2","Normals:0","Normals:1","Normals:2","U:0","U:1","U:2","alpha1","gPc:0","gPc:1","gPc:2","p","pc","gPd:0","gPd:1","gPd:2","U1:0","U1:1","U1:2","p1","pc1","gPd1:0","gPd1:1","gPd1:2","gPc1:0","gPc1:1","gPc1:2","Ux","Uy","Uz","gUx:0","gUx:1","gUx:2","gUy:0","gUy:1","gUy:2","gUz:0","gUz:1","gUz:2","PcxN:0","PcxN:1","PcxN:2","pxN:0","pxN:1","pxN:2","phiU:0","phiU:1","phiU:2","gUzxN","gUxxN","gUyxN","PcxN1:0","PcxN1:1","PcxN1:2","pxN1:0","pxN1:1","pxN1:2","gUzxN1","gUxxN1","Area","Cell Type" ...
% };
%
% %strcmpi(Names,Name)
% returnedDataInd=min(find(strcmpi(Name,Names)));



function  returnedDataInd=    VolI(Names,Name)
	% Names={ ...
	% "t","maxMagU","aAvg","aAvgL","aAvgR","avgUAlpha1_0","avgUAlpha1_1","avgUAlpha1_2","avgUAlpha2_0","avgUAlpha2_1","avgUAlpha2_2","uBack","uFront","Dp","Dpc","pcAvg","ADarcy","S1-alpha","S1-U","S1-vol","S1-f_1","S1-dpddz","S1-dpcdz","S1-dpcdz_1","S1-dpddz_1","S1-viscz","S1-viscz_1","S1-phiu","S1-phiu_1","S1-z1","S1-z2","S2-alpha","S2-U","S2-vol","S2-f_1","S2-dpddz","S2-dpcdz","S2-dpcdz_1","S2-dpddz_1","S2-viscz","S2-viscz_1","S2-phiu","S2-phiu_1","S2-z1","S2-z2","S3-alpha","S3-U","S3-vol","S3-f_1","S3-dpddz","S3-dpcdz","S3-dpcdz_1","S3-dpddz_1","S3-viscz","S3-viscz_1","S3-phiu","S3-phiu_1","S3-z1","S3-z2","S4-alpha","S4-U","S4-vol","S4-f_1","S4-dpddz","S4-dpcdz","S4-dpcdz_1","S4-dpddz_1","S4-viscz","S4-viscz_1","S4-phiu","S4-phiu_1","S4-z1","S4-z2","S5-alpha","S5-U","S5-vol","S5-f_1","S5-dpddz","S5-dpcdz","S5-dpcdz_1","S5-dpddz_1","S5-viscz","S5-viscz_1","S5-phiu","S5-phiu_1","S5-z1","S5-z2" ...
	% };
	returnedDataInd=min(find(strcmp(Name,Names)));

% function  [x_a1,y_a1]= Group(xx,yy,nS_)
	%nS_=20;
	% xx_a1=zeros(1,nS_+1);
	% yy_a1=zeros(1,nS_+1);
	% n_a1=1e-12*ones(1,nS_+1);

	%tmp_mo=0.;
	% for i=1:length(xx)
			% iS_=round(xx(i)*nS_)+1;

			% n_a1(iS_)=n_a1(iS_)+1.;
			% xx_a1(iS_)=xx_a1(iS_)+xx(i);
			% yy_a1(iS_)=yy_a1(iS_)+yy(i);
	% end

	% iS2=1;
	% for iS_=1:nS_
		% if (n_a1(iS_)>0.5)
			% x_a1(iS2)=xx_a1(iS_)/n_a1(iS_);
			% y_a1(iS2)=yy_a1(iS_)/n_a1(iS_);
			% iS2=iS2+1;
		% end
	% end
end




function  xx_a1Avg= groupTimeAvg(xxx,Names,VarName,indexes,delNGroup)
	% VarName
	% index=;
	weight=groupTime(xxx,VolI(Names,["S" num2str(indexes(1)) "-vol"]),delNGroup);
	xx_a1SumWeight=weight;
	xx_a1Sum=weight.*groupTime(xxx,VolI(Names,["S" num2str(indexes(1)) "-" VarName]),delNGroup);
	% mean
	for i=2:length(indexes)
	weight=groupTime(xxx,VolI(Names,["S" num2str(indexes(i)) "-vol"]),delNGroup);
	xx_a1SumWeight=xx_a1SumWeight+weight;
	xx_a1Sum=xx_a1Sum+weight.*groupTime(xxx,VolI(Names,["S" num2str(indexes(i)) "-" VarName]),delNGroup);
	end
	xx_a1Avg=xx_a1Sum./xx_a1SumWeight;
end

function  xx_aMin= groupTimeMin(xxx,Names,VarName,indexes,delNGroup)
	index=VolI(Names,["S" num2str(indexes(1)) "-" VarName]);
	xx_aMin=groupTime(xxx,index,delNGroup);
	for i=2:length(indexes)
		index=VolI(Names,["S" num2str(indexes(i)) "-" VarName]);
		xx_aMin=min(xx_aMin,groupTime(xxx,index,delNGroup));
		minz=min(xx_aMin)
	end
end

function  xx_a1Max= groupTimeMax(xxx,Names,VarName,indexes,delNGroup)
	index=VolI(Names,["S" num2str(indexes(1)) "-" VarName]);
	xx_a1Max=groupTime(xxx,index,delNGroup);
	for i=2:length(indexes)
		index=VolI(Names,["S" num2str(indexes(i)) "-" VarName]);
		xx_a1Max=max(xx_a1Max,groupTime(xxx,index,delNGroup));
		maxz=max(xx_a1Max)
	end
end

function  xx_a1Sum= groupTimeSum(xxx,Names,VarName,indexes,delNGroup)
	index=VolI(Names,["S" num2str(indexes(1)) "-" VarName]);
	xx_a1Sum=groupTime(xxx,index,delNGroup);
	for i=2:length(indexes)
		index=VolI(Names,["S" num2str(indexes(i)) "-" VarName]);
		xx_a1Sum=(xx_a1Sum+groupTime(xxx,index,delNGroup));
		xx_Sum=max(xx_a1Sum)
	end
end

function  xx_a1= groupTime(xxx,index,delNGroup)
	xx=xxx(1:end,index);
	xx_a1=zeros(1,floor(length(xx)/delNGroup));
	n_a1=delNGroup-2;
	if(n_a1>1)
		for i = 1:length(xx_a1)  ,
			ixx2=(i)*delNGroup;
			group=xx(ixx2-delNGroup+1:ixx2);
			xx_a1(i)=(sum(group)-max(group)-min(group))/n_a1 ;
		end
	else
		for i = 1:length(xx_a1)  ,
			ixx2=(i)*delNGroup;
			group=xx(ixx2-delNGroup+1:ixx2);
			xx_a1(i)=mean(group);
		end
	end
	% size(xx_a1)
end




function print2pdflatex (option,file)
	% [fid,name,msg]=mkstemp("./octave-pdflatex-XXXXXX",true)
	% suffix="Berea12_100"
	% drawnow("epscairo",[suffix file ".pdf"],false,[suffix file ".gplt"]);
	% system([" bash bashGnouplotFix " suffix file " " option]);
	% system(["rm  " suffix file ".gplt" ]);
	% system(["convert -density 90 -crop +2+2 " suffix file ".pdf   " suffix file ".png"] );
	% system(sprintf("sed -r -e '2{s|.+|set terminal cairolatex pdf linewidth 4 color|}' %s | gnuplot",name));
	% endfunction
 	drawnow('epscairo',[file '.pdf'],false,[file '.gplt']);
 	system(['convert -density 150 -flatten -crop +2+2 ' file '.pdf ' file '.jpg'] );

end
