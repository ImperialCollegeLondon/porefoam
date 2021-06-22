	Info<< "Reading field pd\n" << endl;
	volScalarField pd(
		IOobject (	"pd", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ),
		mesh );
	
	
	Info<< "Reading field alpha1\n" << endl;
	volScalarField alpha1 (
		IOobject (
			"alpha1", runTime.timeName(), mesh,
			IOobject::MUST_READ, IOobject::NO_WRITE ),
		mesh );
	alpha1.boundaryField().evaluate();
	
//alpha1=fvc::average(fvc::interpolate(alpha1) );
	Info<< "Reading field U\n" << endl;
	volVectorField U (
		IOobject( "U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ),
		mesh );

	#include "./createPhi.H"


	Info<< "Reading transportProperties\n" << endl;
	twoPhaseMixture twoPhaseProperties(U, phi);

	const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
	const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

	const dimensionedScalar& mu1Avg = rho1*average(twoPhaseProperties.nuModel1().nu());
	const dimensionedScalar& mu2Avg = rho2*average(twoPhaseProperties.nuModel2().nu());
	
	Info<< "Laminar viscosities: " <<mu1Avg.value()<<" "<<mu2Avg.value()<< endl;

	// Need to store rho for ddt(rho, U)
	volScalarField rho(
		IOobject("rho", runTime.timeName(), mesh,
				  IOobject::READ_IF_PRESENT ),
		alpha1*rho1 + (scalar(1) - alpha1)*rho2,
		alpha1.boundaryField().types() );
	rho.oldTime();


	// Mass flux
	// Initialisation does not matter because rhoPhi is reset after the
	// alpha1 solution before it is used in the U equation.
	surfaceScalarField rhoPhi(
		IOobject( "rho*phi",  runTime.timeName(),  mesh ),
		rho1*phi );


#ifndef  SINGLE_PHASE
	Info<< "setting up interface \n" << endl;
	// Construct interface from alpha1 distribution
	interfaceProperties interface(alpha1, U, twoPhaseProperties);
	surfaceScalarField deltaS=fvc::snGrad(alpha1);	
#endif
 
	// Construct incompressible turbulence model
	//autoPtr<incompressible::turbulenceModel> turbulence
	//(
		//incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
	//);


//#define errorAnalysis


#ifdef errorAnalysis
	//setRefCell(pd, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

//aq100127
	dimensionedScalar c ( mesh.solutionDict().subDict("errorAnalysis").lookup("c")  );
	volScalarField UAn
	(
		IOobject( "UAn",  runTime.timeName(),  mesh,
			       IOobject::NO_READ,  IOobject::AUTO_WRITE ),
		mesh,
		dimensionedScalar( "zero", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0 ) 
	);


	
	dimensionedScalar h ( mesh.solutionDict().subDict("errorAnalysis").lookup("h")  );
	dimensionedScalar lambdaEA ( mesh.solutionDict().subDict("errorAnalysis").lookup("lambEA")  );
	dimensionedScalar dpPdx ( mesh.solutionDict().subDict("errorAnalysis").lookup("dpPdx")  );
	dimensionedScalar mu ( mesh.solutionDict().subDict("errorAnalysis").lookup("mu")  );
	UAn=(((mesh.C()).component(1)-c)*((mesh.C()).component(1)-c)-h*h-2.*lambdaEA*h)/2/mu*(-dpPdx);
	
	volScalarField Uerror	(
		IOobject( "Uerror",  runTime.timeName(),  mesh,
		          IOobject::NO_READ,  IOobject::AUTO_WRITE		),
	mesh,	0.0	);

#endif 




	label pRefCell = 0;
	scalar pRefValue = 0.0;
	setRefCell(pd, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);



	#include "readGravitationalAcceleration.H"

	surfaceScalarField ghf("ghf", g & mesh.Cf());



	volVectorField UFilter	(
		IOobject( "UFilter",  runTime.timeName(),  mesh	),	
		2.0*(fvc::average(fvc::interpolate(U))),
		"zeroGradient"	);

	label nAlphaCorr( readLabel(mesh.solutionDict().subDict("PIMPLE").lookup("nAlphaCorr")));


	volVectorField gradP=	fvc::grad(pd); 













	IOdictionary meshingDict(
		IOobject( "meshingDict",  runTime.system(),  runTime,
		          IOobject::MUST_READ,  IOobject::NO_WRITE ) );	




	surfaceScalarField muEff
	(
		"muEff",
		twoPhaseProperties.muf()
	);



  #include "createCVs.H"


