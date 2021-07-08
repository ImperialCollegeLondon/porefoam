
//const surfaceVectorField & 	Sf = mesh.Sf();
//const_cast<pointField*>(&mesh.points())->resize(0);
//const_cast<faceList*>(&mesh.faces())->resize(0);

	Info<< "Reading field U\n" << endl;
	volVectorField U
	(	IOobject
		(
			"U",  runTime.timeName(),  mesh,
			IOobject::MUST_READ,  IOobject::NO_WRITE
		),  mesh
	);

	Info<< "Reading field p\n" << endl;
	volScalarField p
	(  IOobject
		(  "p",  runTime.timeName(),  mesh,
			IOobject::MUST_READ,  IOobject::NO_WRITE
		),  mesh
	);





	#include "createPhi.H"

	volScalarField pOldOld=p; 
	volScalarField pOld=p; 


	Info<< "Reading transportProperties\n" << endl;
	IOdictionary transportProperties
	(	IOobject
		( "transportProperties",  runTime.constant(),  runTime,
		IOobject::MUST_READ,  IOobject::NO_WRITE
	)	);

	const dimensionedScalar rho(transportProperties.lookup("rho"));
	const dimensionedScalar nu(transportProperties.lookup("nu"));
	const dimensionedScalar mu(rho*nu);



	scalar cBC =readScalar( pimple.dict().lookup("cBC")  );



	label pRefCell = 0;
	scalar pRefValue = 0.;
	setRefCell(p, pimple.dict(), pRefCell, pRefValue);







	scalar lambda( readScalar(pimple.dict().lookup("lambda")) );



	U.correctBoundaryConditions();



	dimensionedScalar beta ( pimple.dict().lookup("beta")  );

