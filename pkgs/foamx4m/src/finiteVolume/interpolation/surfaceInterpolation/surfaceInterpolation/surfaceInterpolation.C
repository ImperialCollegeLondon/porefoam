/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
	This file is part of foam-extend.

	foam-extend is free software: you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published by the
	Free Software Foundation, either version 3 of the License, or (at your
	option) any later version.

	foam-extend is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "demandDrivenData.H"
#include "coupledFvPatch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(surfaceInterpolation, 0);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::surfaceInterpolation::clearOut()
{
	deleteDemandDrivenData(weightingFactors_);
	deleteDemandDrivenData(deltaCoeffs_);

	orthogonal_ = false;
	deleteDemandDrivenData(correctionVectors_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::surfaceInterpolation::surfaceInterpolation(const fvMesh& fvm)
:
	mesh_(fvm),
	schemesDict_(fvm),
	solutionDict_(fvm),
	weightingFactors_(nullptr),
	deltaCoeffs_(nullptr),
	orthogonal_(false),
	correctionVectors_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::surfaceInterpolation::~surfaceInterpolation()
{
	clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::surfaceScalarField& Foam::surfaceInterpolation::weights() const
{
	if (!weightingFactors_)
	{
		makeWeights();
	}

	return (*weightingFactors_);
}


const Foam::surfaceScalarField& Foam::surfaceInterpolation::deltaCoeffs() const
{
	if (!deltaCoeffs_)
	{
		makeDeltaCoeffs();
	}

	return (*deltaCoeffs_);
}


bool Foam::surfaceInterpolation::orthogonal() const
{
	if (orthogonal_ == false && !correctionVectors_)
	{
		makeCorrectionVectors();
	}

	return orthogonal_;
}


const Foam::surfaceVectorField&
Foam::surfaceInterpolation::correctionVectors() const
{
	if (orthogonal())
	{
		FatalErrorIn("surfaceInterpolation::correctionVectors()")
			<< "cannot return correctionVectors; mesh is orthogonal"
			<< abort(FatalError);
	}

	return (*correctionVectors_);
}


bool Foam::surfaceInterpolation::movePoints()
{
	clearOut();

	return true;
}


void Foam::surfaceInterpolation::makeWeights() const
{
	if (debug)
	{
		Info<< "surfaceInterpolation::makeWeights() : "
			<< "Constructing weighting factors for face interpolation"
			<< endl;
	}

	weightingFactors_ = new surfaceScalarField
	(
		IOobject
		(
			"weightingFactors",
			mesh_.pointsInstance(),
			mesh_
		),
		mesh_,
		dimless
	);
	surfaceScalarField& weightingFactors = *weightingFactors_;


	// Set local references to mesh data
	// (note that we should not use fvMesh sliced fields at this point yet
	//  since this causes a loop when generating weighting factors in
	//  coupledFvPatchField evaluation phase)
	const unallocLabelList& owner = mesh_.owner();
	const unallocLabelList& neighbour = mesh_.neighbour();

	const vectorField& Cf = mesh_.faceCentres();
	const vectorField& C = mesh_.cellCentres();
	const vectorField& Sf = mesh_.faceAreas();

	// ... and reference to the internal field of the weighting factors
	scalarField& w = weightingFactors.internalField();

	forAll (owner, facei)
	{
		// Note: mag in the dot-product.
		// For all valid meshes, the non-orthogonality will be less than
		// 90 deg and the dot-product will be positive.  For invalid
		// meshes (d & s <= 0), this will stabilise the calculation
		// but the result will be poor.
		scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
		scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));
		w[facei] = SfdNei/(SfdOwn + SfdNei);
	}

	forAll (mesh_.boundary(), patchi)
	{
		mesh_.boundary()[patchi].makeWeights
		(
			weightingFactors.boundaryField()[patchi]
		);
	}

	if (debug)
	{
		Info<< "surfaceInterpolation::makeWeights() : "
			<< "Finished constructing weighting factors for face interpolation"
			<< endl;
	}
}


void Foam::surfaceInterpolation::makeDeltaCoeffs() const
{
	if (debug)
	{
		Info<< "surfaceInterpolation::makeDeltaCoeffs() : "
			<< "Constructing differencing factors array for face gradient"
			<< endl;
	}

	// Force the construction of the weighting factors
	// needed to make sure deltaCoeffs are calculated for parallel runs.
	weights();

	deltaCoeffs_ = new surfaceScalarField
	(
		IOobject
		(
			"deltaCoeffs",
			mesh_.pointsInstance(),
			mesh_
		),
		mesh_,
		dimless/dimLength
	);
	surfaceScalarField& DeltaCoeffs = *deltaCoeffs_;


	// Set local references to mesh data
	const volVectorField& C = mesh_.C();
	const unallocLabelList& owner = mesh_.owner();
	const unallocLabelList& neighbour = mesh_.neighbour();
	const surfaceVectorField& Sf = mesh_.Sf();
	const surfaceScalarField& magSf = mesh_.magSf();

	forAll (owner, facei)
	{
		vector delta = C[neighbour[facei]] - C[owner[facei]];
		vector unitArea = Sf[facei]/magSf[facei];

		// Standard cell-centre distance form
		//DeltaCoeffs[facei] = (unitArea & delta)/magSqr(delta);

		// Slightly under-relaxed form
		//DeltaCoeffs[facei] = 1.0/mag(delta);

		// More under-relaxed form
		//DeltaCoeffs[facei] = 1.0/(mag(unitArea & delta) + VSMALL);

		// Stabilised form for bad meshes
		DeltaCoeffs[facei] = 1.0/max(unitArea & delta, 0.05*mag(delta));
	}

	forAll (DeltaCoeffs.boundaryField(), patchi)
	{
		mesh_.boundary()[patchi].makeDeltaCoeffs
		(
			DeltaCoeffs.boundaryField()[patchi]
		);
	}
}


void Foam::surfaceInterpolation::makeCorrectionVectors() const
{
	if (debug)
	{
		Info<< "surfaceInterpolation::makeCorrectionVectors() : "
			<< "Constructing non-orthogonal correction vectors"
			<< endl;
	}

	correctionVectors_ = new surfaceVectorField
	(
		IOobject
		(
			"correctionVectors",
			mesh_.pointsInstance(),
			mesh_
		),
		mesh_,
		dimless
	);
	surfaceVectorField& corrVecs = *correctionVectors_;

	// Set local references to mesh data
	const volVectorField& C = mesh_.C();
	const unallocLabelList& owner = mesh_.owner();
	const unallocLabelList& neighbour = mesh_.neighbour();
	const surfaceVectorField& Sf = mesh_.Sf();
	const surfaceScalarField& magSf = mesh_.magSf();
	const surfaceScalarField& DeltaCoeffs = deltaCoeffs();

	forAll (owner, facei)
	{
		vector unitArea = Sf[facei]/magSf[facei];
		vector delta = C[neighbour[facei]] - C[owner[facei]];

		// If non-orthogonality is over 90 deg, kill correction vector
		// HJ, 27/Feb/2011
		corrVecs[facei] = pos(unitArea & delta)*
			(unitArea - delta*DeltaCoeffs[facei]);
	}

	// Boundary correction vectors set to zero for boundary patches
	// and calculated consistently with internal corrections for
	// coupled patches

	forAll (corrVecs.boundaryField(), patchI)
	{
		mesh_.boundary()[patchI].makeCorrVecs
		(
			corrVecs.boundaryField()[patchI]
		);
	}

	scalar MaxNonOrthog = 0.0;

	// Calculate the non-orthogonality for meshes with 1 face or more
	if (returnReduce(magSf.size(), sumOp<label>()) > 0)
	{
		MaxNonOrthog =
			asin
			(
				min
				(
					max(mag(corrVecs)).value(),
					1.0
				)
			)*180.0/mathematicalConstant::pi;
	}

	if (debug)
	{
		Info<< "surfaceInterpolation::makeCorrectionVectors() : "
			<< "maximum non-orthogonality = " << MaxNonOrthog << " deg."
			<< endl;
	}

	if (MaxNonOrthog < 5)
	{
		orthogonal_ = true;
		deleteDemandDrivenData(correctionVectors_);
	}
	else
	{
		orthogonal_ = false;
	}

	if (debug)
	{
		Info<< "surfaceInterpolation::makeCorrectionVectors() : "
			<< "Finished constructing non-orthogonal correction vectors"
			<< endl;
	}
}


// ************************************************************************* //
