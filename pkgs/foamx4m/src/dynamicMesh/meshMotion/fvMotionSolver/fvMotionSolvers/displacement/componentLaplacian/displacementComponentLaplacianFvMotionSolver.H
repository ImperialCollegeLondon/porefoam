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

Class
	Foam::displacementComponentLaplacianFvMotionSolver

Description
	Mesh motion solver for an fvMesh.  Based on solving the cell-centre
	Laplacian for the given component of the motion displacement.

SourceFiles
	displacementComponentLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementComponentLaplacianFvMotionSolver_H
#define displacementComponentLaplacianFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;


class displacementComponentLaplacianFvMotionSolver
:
	public fvMotionSolver
{
	// Private data

		//- The component name to solve for
		word cmptName_;

		//- The component to solve for
		direction cmpt_;

		//- Reference point field for this component
		scalarField points0_;

		//- Point motion field
		mutable pointScalarField pointDisplacement_;

		//- Cell-centre motion field
		mutable volScalarField cellDisplacement_;

		//- Optionally read point-position field. Used only for position
		//  boundary conditions.
		mutable autoPtr<pointVectorField> pointLocation_;

		//- Diffusivity used to control the motion
		autoPtr<motionDiffusivity> diffusivityPtr_;

		//- Frozen points (that are not on patches). -1 or points that are
		//  fixed to be at points0_ location
		label frozenPointsZone_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		displacementComponentLaplacianFvMotionSolver
		(
			const displacementComponentLaplacianFvMotionSolver&
		);

		//- Disallow default bitwise assignment
		void operator=(const displacementComponentLaplacianFvMotionSolver&);

		//- Return the component corresponding to the given component name
		direction cmpt(const word& cmptName) const;


public:

	//- Runtime type information
	TypeName("displacementComponentLaplacian");


	// Constructors

		//- Construct from polyMesh and data stream
		displacementComponentLaplacianFvMotionSolver
		(
			const polyMesh&,
			Istream& msData
		);


	// Destructor

		~displacementComponentLaplacianFvMotionSolver();


	// Member Functions

		//- Return reference to the reference field
		const scalarField& points0() const
		{
			return points0_;
		}

		//- Non-const access to the cellDisplacement in order to allow
		//  changes to the boundary motion
		volScalarField& cellDisplacement()
		{
			return cellDisplacement_;
		}

		//- Return point location obtained from the current motion field
		virtual tmp<pointField> curPoints() const;

		//- Solve for motion
		virtual void solve();

		//- Update topology
		virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
