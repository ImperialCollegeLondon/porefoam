/*-------------------------------------------------------------------------*\
This code is part of poreFOAM, a suite of codes written using OpenFOAM
for direct simulation of flow at the pore scale. 	
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.



The code has been developed by Ali Qaseminejad Raeini as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
 
 Description:
	post-processing code Calculates energy losses for control 
	volumes given in postprossingDict file in $case/system directory. 
	
	Note: Parallel case post-processing is added to the code, no need to 
	reconstract parallel cases anymore.
\*-------------------------------------------------------------------------*/






//int PostProcess(const fvMesh& mesh)
if(CVBounds1.size())
{
	Info<<endl;


	snapShot_t data(CVBounds1.size());

	data.t=runTime.value();

	volScalarField Uz=U.component(0);
	volVectorField gUz=fvc::grad(Uz);
	volScalarField gradPz=gradP.component(0);
	volScalarField gPcz=gPc.component(0);


	volScalarField ViscERate=(fvc::laplacian(muEff,U) & U);
	volScalarField gPcERate=gPc & U;
	volScalarField gradPERate=-gradP & U;
	volScalarField phiERate=-fvc::div(rhoPhi, U) & U ;



	scalar cPc2=max(0.5*(1.0-interface.pcThicknessFactor()),.45);
	volScalarField alpha1Sh=(1/(cPc2))*(max(alpha1,(1-cPc2))-(1-cPc2));
	volScalarField alpha2Sh=1-(1/(cPc2))*min(alpha1,(cPc2));
	volScalarField alphaSh=(1./(1-2*0.20))*(curtailBADOFSET(alpha1,0.20)-0.20);



	data.aAvg=average(alpha1).value();
	data.aAvgL=gSum(alpha1.boundaryField()[iBack]*alpha1.boundaryField()[iBack].patch().magSf())
				/(gSum(alpha1.boundaryField()[iBack].patch().magSf())+1.0e-24);
	data.aAvgR=gSum(alpha1.boundaryField()[iFront]*alpha1.boundaryField()[iFront].patch().magSf())
				/(gSum(alpha1.boundaryField()[iFront].patch().magSf())+1.0e-24);


	data.pcAvg = (sum(alpha1Sh*pc)/(mag(sum(alpha1Sh))+1e-12)).value();
	data.pcAvg -= (sum(alpha2Sh*pc)/(mag(sum(alpha2Sh))+1e-12)).value();



	data.maxMagU= max(mag(U)).value();
	data.avgUAlpha1=(sum(alpha1Sh*U)/(mag(sum(alpha1Sh))+1e-12)).value();
	data.avgUAlpha2=(sum(alpha2Sh*U)/(mag(sum(alpha2Sh))+1e-12)).value();
	Info <<"___________  Sw:"<<1.0-data.aAvg <<"   Pc:"<<data.pcAvg <<"   U1:"<<data.avgUAlpha1<<"   U2:"<<data.avgUAlpha2<<"   ____________"<< endl;

	data.QIn=gSum(U.boundaryField()[iBack]&U.boundaryField()[iBack].patch().Sf());
	data.QOut=gSum(U.boundaryField()[iFront]&U.boundaryField()[iFront].patch().Sf());
	
	
	data.Dp=(gAverage(pd.boundaryField()[iBack])-gAverage(pd.boundaryField()[iFront]));
	data.Dpc=(gAverage(pc.boundaryField()[iBack])-gAverage(pc.boundaryField()[iFront]));
	Info <<"DPc"<< data.Dpc << endl;


	data.ADarcy=(max(mesh.C().component(2)).value()-min(mesh.C().component(2)).value())*(max(mesh.C().component(1)).value()-min(mesh.C().component(1)).value());




	const vectorField & C =	mesh.C().internalField();

	for (int iSam=0;iSam<regionVoxelValues.size();iSam++)   
	{ 
		snapShot_zt & slice = data.slices[iSam];
		double xSam1=CVBounds1[iSam];
		double xSam2=CVBounds2[iSam];

		Info<<iSam<< "   Del x :"<< xSam1 <<" to "<< xSam2<< "   ";		

		//clip1=0.0;
		clip=0.0;
		forAll(Cx,c) 
		{ 
				if (PPRegions[c]==regionVoxelValues[iSam]&& 
					(C[c][0]>=CVBounds1[iSam] && C[c][0]<=CVBounds2[iSam] && 
					 C[c][1]>=CVyBounds1[iSam] && C[c][1]<=CVyBounds2[iSam] && 
					 C[c][2]>=CVzBounds1[iSam] && C[c][2]<=CVzBounds2[iSam] ))
				   clip[c]=1.0;
		}




		slice.mu1()=mu1Avg.value();
		slice.mu2()=mu2Avg.value();

		slice.x1()=xSam1;		
		slice.x2()=xSam2;
		scalar delX=xSam2-xSam1;
		Info<<"   ->         "<<slice.x1()<<"  "<<slice.x2()<<endl;		


		slice.zDropAvg()=gAverage((Cx.internalField())*(clip*(1-alphaSh))->internalField());
		slice.zDrop1()=xSam2-gMax((xSam2-Cx.internalField())*(clip*(1-alphaSh))->internalField());
		slice.zDrop2()=gMax((Cx*clip*(1-alphaSh))->internalField());


		volScalarField clipAlpha=alphaSh*clip;
		volScalarField clipVC=clip*VC;
		const fvBoundaryMesh& boundary = mesh.boundary();

		forAll(boundary, patchi)
		{
			clip.boundaryField()[patchi]=0.0;
			clip.boundaryField()[patchi]==0.0;
			clipAlpha.boundaryField()[patchi]=0.0;
			clipAlpha.boundaryField()[patchi]==0.0;
		}
		surfaceScalarField BSurafce=fvcDelta(clip);
		surfaceScalarField AlphaInternalBoundary=fvcDelta(clipAlpha);
		volScalarField clipAlphaVC=clipAlpha*VC;



		scalar vol=sum(clipVC).value();

		
		if(mag(vol)<1.0e-31) 
		{
			Info<<"Error in post-processing: CV volume is zero: label "<<regionVoxelValues[iSam]<<endl;
			vol=1.0e-31;
		}
		slice.vol()=vol;

		slice.alpha()=sum(clipAlphaVC).value()/vol;
		slice.q()=sum(Uz *clipVC).value()/vol; if(mag(slice.q())<1.0e-32) slice.q()=1.0e-32;
		slice.f_1()=sum(Uz*clipAlphaVC).value()/vol/slice.q();
		slice.dpddz()=-sum(gradPz *clipVC).value()/vol;
		slice.dpcdz()=sum(gPcz *clipVC).value()/vol;
		slice.dpddz_1()=-sum(gradPz*clipAlphaVC).value()/vol;
		slice.dpcdz_1()=sum(gPcz*clipAlphaVC).value()/vol;

		#include  "calcdEdts.H"


		slice.viscz()=viscF[2]/vol;
		slice.viscz_1()=viscF1[2]/vol; 
		slice.phiu()=rhoPhiU[2]/vol;
		slice.phiu_1()=rhoPhiU1[2]/vol;
		slice.delPdelZ()=gPd[2]/vol;



		slice.viscInterf_1()=viscInterf_1[2]/vol;

		slice.viscE()=sum( ViscERate  *clipVC).value()/vol;
		slice.viscE_1()=sum( ViscERate * clipAlphaVC).value()/vol;
		slice.dpEd()=sum( gradPERate *clipVC).value()/vol;
		slice.dpEd_1()=sum(gradPERate * clipAlphaVC).value()/vol;
		slice.dpEc()=sum( gPcERate *clipVC).value()/vol;
		slice.dpEc_1()=sum( gPcERate * clipAlphaVC).value()/vol;
		slice.phiE()=sum( phiERate *clipVC).value()/vol;
		slice.phiE_1()=sum( phiERate * clipAlphaVC).value()/vol;


		Info<< "viscE phiE dpEd dpEc:  "
			 << slice.viscE()*vol/(slice.q()*vol/delX)<< "    " 
			 << slice.phiE()*vol/(slice.q()*vol/delX)   << "      " 
			 << slice.dpEd()*vol/(slice.q()*vol/delX)   << "    " 
			 << slice.dpEc()*vol/(slice.q()*vol/delX)   <<endl;


		slice.dpddz_1()=-sum(gradPz*clipAlphaVC).value()/vol;
		slice.dpcdz_1()=sum(gPcz*clipAlphaVC).value()/vol;

		slice.Pc()=sum(pc*alpha1Sh *clipVC).value()/max(sum(alpha1Sh *clipVC).value(),5.0e-32);
		slice.Pc()-=sum(pc*alpha2Sh *clipVC).value()/max(sum(alpha2Sh *clipVC).value(),5.0e-32);


	}

	if (Pstream::master())
	{
		data.write(data_out_for_plot);
		Info<<"wrote data_out_for_plot"<<endl;
	}

	Info<<endl;
} 



