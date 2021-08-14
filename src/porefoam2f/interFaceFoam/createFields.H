
	Info<< "Reading field pd\n" << endl;
	volScalarField pd(
		IOobject( "pd", runTime.timeName(), mesh,
			       IOobject::MUST_READ, IOobject::AUTO_WRITE ), 
		mesh );


	Info<< "Reading field alpha1\n" << endl;
	volScalarField alpha1 (
		IOobject("alpha1", runTime.timeName(), mesh,
				IOobject::MUST_READ, IOobject::AUTO_WRITE ),
		mesh );
	alpha1=curtailBADOFSET(alpha1,0.);
	alpha1.boundaryFieldRef().evaluate();


	Info<< "Reading field U\n" << endl;
	volVectorField U (
		IOobject( "U",  runTime.timeName(),  mesh,
			       IOobject::MUST_READ,  IOobject::AUTO_WRITE ),
		mesh );

	#include "createPhi.H"


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



	Info<< "setting up interface \n" << endl;
	// Construct interface from alpha1 distribution
	interfaceProps interface(alpha1, U, twoPhaseProperties);

 
	// Construct incompressible turbulence model
	//autoPtr<incompressible::turbulenceModel> turbulence
	//(
		//incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
	//);





	label pRefCell, pRefCellOrig = 0;
	scalar pRefValue, pRefValueOrig = 0.;
	setRefCell(pd, mesh.solutionDict().subDict("PIMPLE"), pRefCellOrig, pRefValueOrig);
#define SET_REF_(rfc, rfv, sfx, alpha) \
	for(rfc = rfc##sfx; rfc<alpha.size()-1; ++rfc) if (alpha[rfc]<0.01)  break;   \
	rfv=0.99*pd[rfc]+0.01*rfv##sfx

	SET_REF_(pRefCell,pRefValue,Orig,alpha1);

	#include "readGravitationalAcceleration.H"

	surfaceScalarField ghf("ghf", g & mesh.Cf());



	volVectorField UFilter	(
		IOobject( "UFilter",  runTime.timeName(),  mesh	),	
		2.*(fvc::average(fvc::interpolate(U))),
		"zeroGradient"	);

	label nAlphaCorr( readLabel(mesh.solutionDict().subDict("PIMPLE").lookup("nAlphaCorr")));


	const surfaceScalarField& Internalfaces1 = interface.Internalfaces1();

	volScalarField pFilter
	(
		IOobject( "pFilter",  runTime.timeName(),  mesh	),	
		2.*(fvc::average(fvc::interpolate(pd))),	"zeroGradient"
	);



	dimensionedScalar deltaN(1e-8/Foam::cbrt(gAverage(mesh.V())));
	Info<<"<<deltaN.value()="<<deltaN.value()<<endl;

	//pd.correctBoundaryConditions();
	U.correctBoundaryConditions();



	dimensionedScalar beta ( mesh.solutionDict().subDict("PIMPLE").lookup("beta")  );


	scalar cBC =readScalar( mesh.solutionDict().subDict("PIMPLE").lookup("cBC")  );

	scalar nSlip(readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("nSlip")));
	scalar cSSlip(readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("cSSlip")));

	scalar lambdaS(readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("lambdaS")));
	scalar lambda(readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("lambda")));



	scalar deltaTFact=1.;
	scalar CoNumOld=0.;

	label filterU (  mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("filterU", label(1))  );
	label filterPd (  mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("filterPd", label(1))  );
	scalar UBoundByUSmoothFactor 	(  mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("UBoundByUSmoothFactor", scalar(2.))  );
	scalar pBoundBypSmoothFactor	(  mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("pBoundBypSmoothFactor", scalar(1.3))  );

	surfaceScalarField sgPc
	(	IOobject( "sgPcRelaxed",  runTime.timeName(),  mesh,
			       IOobject::NO_READ,  IOobject::AUTO_WRITE ),
		interface.sgPc()
	);



	volVectorField gradP=	fvc::grad(pd); 









	IOdictionary meshingDict(
		IOobject( "meshingDict",  runTime.system(),  runTime,
		          IOobject::MUST_READ,  IOobject::NO_WRITE ) );	




	surfaceScalarField muEff
	(
		"muEff",
		twoPhaseProperties.muf()
	);





