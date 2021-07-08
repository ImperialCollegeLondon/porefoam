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
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@imperial.ac.uk

 Description:
	post-processing code Calculates energy losses for control 
	volumes given in postprossingDict file in $case/system directory. 
	
\*-------------------------------------------------------------------------*/






//int PostProcess(const fvMesh& mesh)
{
	Info<<endl;



	surfaceScalarField muEff(IOobject ("muEff",  runTime.timeName(), mesh), mesh, rho*nu);
	#include  "correctmuEff.H"







	//int counter = 0;
	vector uTotal(0.,0.,0.);

//forAll(timeDirs, timeI)
//if (timeDirs[timeI]==0.0671834){
//runTime.setTime(timeDirs[timeI], timeI);

	snapShot_t data(CVBounds1.size());

	data.t=runTime.value();

	volVectorField gradP=	fvc::grad(p); 

	volScalarField Uz=U.component(0);
	volVectorField gUz=fvc::grad(Uz);
	volScalarField gradPz=gradP.component(0);


	volScalarField ViscERate=(fvc::laplacian(muEff,U) & U);
	volScalarField gradPERate=-gradP & U;
	volScalarField phiERate=-fvc::div(rho*phi, U) & U ;







	data.maxMagU= max(mag(U)).value();

	data.QIn=gSum(U.boundaryField()[iBack]&U.boundaryField()[iBack].patch().Sf());
	data.QOut=gSum(U.boundaryField()[iFront]&U.boundaryField()[iFront].patch().Sf());
	
	
	data.Dp=(gAverage(p.boundaryField()[iBack])-gAverage(p.boundaryField()[iFront]));

	const volVectorField & C =	mesh.C();

	data.ADarcy=(max(C.component(2)).value()-min(C.component(2)).value())*(max(C.component(1)).value()-min(C.component(1)).value());





	for (int iSam=0;iSam<CVxlVals.size();iSam++)   
	{
		snapShot_zt & slice = data.slices[iSam];
		double xSam1=CVBounds1[iSam];
		double xSam2=CVBounds2[iSam];

		Info<<iSam<< "   Del x :"<< xSam1 <<" to "<< xSam2<< "   ";		

		volScalarField clip
		(	IOobject( "clip", runTime.timeName(), mesh),
			mesh,	dimensionedScalar("clip",dimless,0.),	"fixedValue"
		);
		forAll(clip,c) 
		{ 
				if (PPRegions[c]==CVxlVals[iSam]&& 
					(C[c][0]>=CVBounds1[iSam] && C[c][0]<=CVBounds2[iSam] && 
					 C[c][1]>=CVyBounds1[iSam] && C[c][1]<=CVyBounds2[iSam] && 
					 C[c][2]>=CVzBounds1[iSam] && C[c][2]<=CVzBounds2[iSam] ))
					clip[c]=1.;
		}


		slice.x1()=xSam1;		
		slice.x2()=xSam2;
		scalar delX=xSam2-xSam1;
		Info<<"   ->         "<<slice.x1()<<"  "<<slice.x2()<<endl;		




		volScalarField clipVC=clip*VC;
		const fvBoundaryMesh& boundary = mesh.boundary();

		forAll(boundary, bi)
		{
			clip.boundaryField()[bi]=0.;
			clip.boundaryField()[bi]==0.;
		}



		scalar vol=sum(clipVC).value();
		
		if(mag(vol)<1e-31) 
		{
			Info<<"Error in post-processing: CV volume is zero: label "<<CVxlVals[iSam]<<endl;
			vol=1e-31;
		}
		slice.vol()=vol;

		slice.q()=sum(Uz *clipVC).value()/vol; if(mag(slice.q())<1e-32) slice.q()=1e-32;
		slice.dpddz()=-sum(gradPz *clipVC).value()/vol;

		//#include  "calcdEdts.H"
		vector viscE(0.,0.,0.);
		vector gPdE(0.,0.,0.);
		vector rhoPhiUE(0.,0.,0.);
		//{

		surfaceScalarField BSurafce=fvcDelta(clip);

		volVectorField::Boundary& Ubf = U.boundaryField();

		forAll(boundary, bi) if (isA<fixedValueFvPatchVectorField>(Ubf[bi]))
		{

			scalarField muEffbp= muEff.boundaryField()[bi];

			scalarField clipbp=    BSurafce.boundaryField()[bi];

			vectorField bGradUp=Ubf[bi].snGrad();//(Ubf[bi].patchInternalField()-Ubf[bi])*boundary[bi].deltaCoeffs();
			bGradUp=bGradUp-boundary[bi].Sf()*(bGradUp&boundary[bi].Sf())/boundary[bi].magSf()/boundary[bi].magSf();
			viscE-=gSum(          clipbp*muEffbp*bGradUp	*boundary[bi].magSf() );


			gPdE+=gSum(clipbp*p.boundaryField()[bi]	*boundary[bi].Sf() );

			rhoPhiUE +=gSum( (clipbp*phi.boundaryField()[bi]*fvc::interpolate(rho*U)->boundaryField()[bi]) );
		}

		surfaceVectorField gIBC= fvc::snGrad(U);
		gIBC=gIBC-mesh.Sf()*(gIBC&mesh.Sf())/mesh.magSf()/mesh.magSf();

		viscE-=gSum( (BSurafce*muEff*gIBC	*mesh.magSf())->internalField() );

		gPdE+=gSum( ( BSurafce*fvc::interpolate(p)	*mesh.Sf())->internalField() );	

		rhoPhiUE +=gSum( (BSurafce*phi*fvc::interpolate(rho*U))->internalField() );
		//}



		slice.viscz()=viscE[2]/vol;
		slice.phiu()=rhoPhiUE[2]/vol;
		slice.delPdelX()=gPdE[2]/vol;




		slice.viscE()=sum( ViscERate  *clipVC).value()/vol;
		slice.dpdE()=sum( gradPERate *clipVC).value()/vol;
		slice.phiE()=sum( phiERate *clipVC).value()/vol;


		Info<< "viscE phiE dpdE :  "
			 << slice.viscE()*vol/(slice.q()*vol/delX)<< "    " 
			 << slice.phiE()*vol/(slice.q()*vol/delX)   << "      " 
			 << slice.dpdE()*vol/(slice.q()*vol/delX)   <<endl;



	}

	if (Pstream::master())
	{
		data.write(data_out_for_plot);
		Info<<"wrote grads_SP.csv"<<endl;
	}
//}

	Info<<endl;
} 




