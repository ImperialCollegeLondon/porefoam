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

#include "velocityComponentLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(velocityComponentLaplacianFvMotionSolver, 0);

	addToRunTimeSelectionTable
	(
		fvMotionSolver,
		velocityComponentLaplacianFvMotionSolver,
		dictionary
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityComponentLaplacianFvMotionSolver::
velocityComponentLaplacianFvMotionSolver
(
	const polyMesh& mesh,
	Istream& msData
)
:
	fvMotionSolver(mesh),
	cmptName_(msData),
	cmpt_(0),
	pointMotionU_
	(
		IOobject
		(
			"pointMotionU" + cmptName_,
			fvMesh_.time().timeName(),
			fvMesh_,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		pointMesh::New(fvMesh_)
	),
	cellMotionU_
	(
		IOobject
		(
			"cellMotionU" + cmptName_,
			mesh.time().timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		fvMesh_,
		dimensionedScalar
		(
			"cellMotionU",
			pointMotionU_.dimensions(),
			0
		),
		cellMotionBoundaryTypes<scalar>(pointMotionU_.boundaryField())
	),
	diffusivityPtr_
	(
		motionDiffusivity::New(*this, lookup("diffusivity"))
	)
{
	if (cmptName_ == "x")
	{
		cmpt_ = vector::X;
	}
	else if (cmptName_ == "y")
	{
		cmpt_ = vector::Y;
	}
	else if (cmptName_ == "z")
	{
		cmpt_ = vector::Z;
	}
	else
	{
		FatalErrorIn
		(
			"velocityComponentLaplacianFvMotionSolver::"
			"velocityComponentLaplacianFvMotionSolver"
			"(const polyMesh& mesh, Istream& msData)"
		)   << "Given component name " << cmptName_ << " should be x, y or z"
			<< exit(FatalError);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityComponentLaplacianFvMotionSolver::
~velocityComponentLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityComponentLaplacianFvMotionSolver::curPoints() const
{
	volPointInterpolation::New(fvMesh_).interpolate
	(
		cellMotionU_,
		pointMotionU_
	);

	tmp<pointField> tcurPoints(new pointField(fvMesh_.points()));

	tcurPoints().replace
	(
		cmpt_,
		tcurPoints().component(cmpt_)
	  + fvMesh_.time().deltaT().value()*pointMotionU_.internalField()
	);

	twoDCorrectPoints(tcurPoints());

	return tcurPoints;
}


void Foam::velocityComponentLaplacianFvMotionSolver::solve()
{
	// The points have moved so before interpolation update
	// the fvMotionSolver accordingly
	movePoints(fvMesh_.points());

	diffusivityPtr_->correct();
	pointMotionU_.boundaryField().updateCoeffs();

	Foam::solve
	(
		fvm::laplacian
		(
			diffusivityPtr_->operator()(),
			cellMotionU_,
			"laplacian(diffusivity,cellMotionU)"
		)
	);
}


void Foam::velocityComponentLaplacianFvMotionSolver::updateMesh
(
	const mapPolyMesh& mpm
)
{
	fvMotionSolver::updateMesh(mpm);

	// Update diffusivity. Note two stage to make sure old one is de-registered
	// before creating/registering new one.
	diffusivityPtr_.reset(nullptr);
	diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
