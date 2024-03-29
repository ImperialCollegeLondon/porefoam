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

#include "displacementLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(displacementLaplacianFvMotionSolver, 0);

	addToRunTimeSelectionTable
	(
		fvMotionSolver,
		displacementLaplacianFvMotionSolver,
		dictionary
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLaplacianFvMotionSolver::displacementLaplacianFvMotionSolver
(
	const polyMesh& mesh,
	Istream& is
)
:
	displacementFvMotionSolver(mesh, is),
	pointDisplacement_
	(
		IOobject
		(
			"pointDisplacement",
			fvMesh_.time().timeName(),
			fvMesh_,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		pointMesh::New(fvMesh_)
	),
	cellDisplacement_
	(
		IOobject
		(
			"cellDisplacement",
			mesh.time().timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		fvMesh_,
		dimensionedVector
		(
			"cellDisplacement",
			pointDisplacement_.dimensions(),
			vector::zero
		),
		cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
	),
	pointLocation_(nullptr),
	diffusivityPtr_
	(
		motionDiffusivity::New(*this, lookup("diffusivity"))
	),
	frozenPointsZone_
	(
		found("frozenPointsZone")
	  ? fvMesh_.pointZones().findZoneID(lookup("frozenPointsZone"))
	  : -1
	)
{
	IOobject io
	(
		"pointLocation",
		fvMesh_.time().timeName(),
		fvMesh_,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	);

	if (debug)
	{
		Info<< "displacementLaplacianFvMotionSolver:" << nl
			<< "    diffusivity	   : " << diffusivityPtr_().type() << nl
			<< "    frozenPoints zone : " << frozenPointsZone_ << endl;
	}


	if (io.headerOk())
	{
		pointLocation_.reset
		(
			new pointVectorField
			(
				io,
				pointMesh::New(fvMesh_)
			)
		);

		if (debug)
		{
			Info<< "displacementLaplacianFvMotionSolver :"
				<< " Read pointVectorField "
				<< io.name() << " to be used for boundary conditions on points."
				<< nl
				<< "Boundary conditions:"
				<< pointLocation_().boundaryField().types() << endl;
		}
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementLaplacianFvMotionSolver::
~displacementLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementLaplacianFvMotionSolver::curPoints() const
{
	volPointInterpolation::New(fvMesh_).interpolate
	(
		cellDisplacement_,
		pointDisplacement_
	);

	if (pointLocation_.valid())
	{
		if (debug)
		{
			Info<< "displacementLaplacianFvMotionSolver : applying "
				<< " boundary conditions on " << pointLocation_().name()
				<< " to new point location."
				<< endl;
		}

		pointLocation_().internalField() =
			points0()
		  + pointDisplacement_.internalField();

		pointLocation_().correctBoundaryConditions();

		// Implement frozen points
		if (frozenPointsZone_ != -1)
		{
			const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

			forAll(pz, i)
			{
				pointLocation_()[pz[i]] = points0()[pz[i]];
			}
		}

		twoDCorrectPoints(pointLocation_().internalField());

		return tmp<pointField>(pointLocation_().internalField());
	}
	else
	{
		tmp<pointField> tcurPoints
		(
			points0() + pointDisplacement_.internalField()
		);

		// Implement frozen points
		if (frozenPointsZone_ != -1)
		{
			const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

			forAll(pz, i)
			{
				tcurPoints()[pz[i]] = points0()[pz[i]];
			}
		}

		twoDCorrectPoints(tcurPoints());

		return tcurPoints;
	}
}


void Foam::displacementLaplacianFvMotionSolver::solve()
{
	// The points have moved so before interpolation update
	// the fvMotionSolver accordingly
	movePoints(fvMesh_.points());

	diffusivityPtr_->correct();
	pointDisplacement_.boundaryField().updateCoeffs();

	Foam::solve
	(
		fvm::laplacian
		(
			diffusivityPtr_->operator()(),
			cellDisplacement_,
			"laplacian(diffusivity,cellDisplacement)"
		)
	);
}


void Foam::displacementLaplacianFvMotionSolver::updateMesh
(
	const mapPolyMesh& mpm
)
{
	displacementFvMotionSolver::updateMesh(mpm);

	// Update diffusivity. Note two stage to make sure old one is de-registered
	// before creating/registering new one.
	diffusivityPtr_.reset(nullptr);
	diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
