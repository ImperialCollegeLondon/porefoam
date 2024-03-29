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
	Foam::velocityComponentLaplacianFvMotionSolver

Description
	Mesh motion solver for an fvMesh.  Based on solving the cell-centre
	Laplacian for the given component of the motion velocity.

SourceFiles
	velocityComponentLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityComponentLaplacianFvMotionSolver_H
#define velocityComponentLaplacianFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;


class velocityComponentLaplacianFvMotionSolver
:
	public fvMotionSolver
{
	// Private data

		//- The component name to solve for
		word cmptName_;

		//- The component to solve for
		direction cmpt_;

		//- Point motion field
		mutable pointScalarField pointMotionU_;

		//- Cell-centre motion field
		mutable volScalarField cellMotionU_;

		//- Diffusivity used to control the motion
		autoPtr<motionDiffusivity> diffusivityPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		velocityComponentLaplacianFvMotionSolver
		(
			const velocityComponentLaplacianFvMotionSolver&
		);

		//- Disallow default bitwise assignment
		void operator=(const velocityComponentLaplacianFvMotionSolver&);


public:

	//- Runtime type information
	TypeName("velocityComponentLaplacian");


	// Constructors

		//- Construct from polyMesh and data stream (provides component)
		velocityComponentLaplacianFvMotionSolver
		(
			const polyMesh&,
			Istream& msData
		);


	// Destructor

		virtual ~velocityComponentLaplacianFvMotionSolver();


	// Member Functions

		//- Non-const access to the cellMotionU in order to allow changes
		//  to the boundary motion
		volScalarField& cellMotionU()
		{
			return cellMotionU_;
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
