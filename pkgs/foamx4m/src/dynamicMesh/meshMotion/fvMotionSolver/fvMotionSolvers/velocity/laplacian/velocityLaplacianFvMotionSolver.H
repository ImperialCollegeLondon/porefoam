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
	Foam::velocityLaplacianFvMotionSolver

Description
	Mesh motion solver for an fvMesh.  Based on solving the cell-centre
	Laplacian for the motion velocity.

SourceFiles
	velocityLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityLaplacianFvMotionSolver_H
#define velocityLaplacianFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;


class velocityLaplacianFvMotionSolver
:
	public fvMotionSolver
{
	// Private data

		//- Point motion field
		mutable pointVectorField pointMotionU_;

		//- Cell-centre motion field
		mutable volVectorField cellMotionU_;

		//- Diffusivity used to control the motion
		autoPtr<motionDiffusivity> diffusivityPtr_;

		//- Use leastSquares volToPoint interpolaltion
		Switch leastSquaresVolPoint_;

	// Private Member Functions

		//- Disallow default bitwise copy construct
		velocityLaplacianFvMotionSolver
		(
			const velocityLaplacianFvMotionSolver&
		);

		//- Disallow default bitwise assignment
		void operator=(const velocityLaplacianFvMotionSolver&);


public:

	//- Runtime type information
	TypeName("velocityLaplacian");


	// Constructors

		//- Construct from polyMesh and data stream
		velocityLaplacianFvMotionSolver
		(
			const polyMesh&,
			Istream& msDataUnused
		);


	// Destructor

		virtual ~velocityLaplacianFvMotionSolver();


	// Member Functions

		//- Return reference to the point motion velocity field
		pointVectorField& pointMotionU()
		{
			return pointMotionU_;
		}

		//- Return const reference to the point motion velocity field
		const pointVectorField& pointMotionU() const
		{
			return pointMotionU_;
		}

		//- Return reference to the cell motion velocity field
		volVectorField& cellMotionU()
		{
			return cellMotionU_;
		}

		//- Return const reference to the cell motion velocity field
		const volVectorField& cellMotionU() const
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
