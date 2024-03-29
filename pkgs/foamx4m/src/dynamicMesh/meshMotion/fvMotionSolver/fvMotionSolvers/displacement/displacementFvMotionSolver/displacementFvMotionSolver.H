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
	Foam::displacementFvMotionSolver.H

Description
	Base class for fvMotionSolvers which calculate displacement.

SourceFiles
	displacementFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementFvMotionSolver_H
#define displacementFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class displacementFvMotionSolver
:
	public fvMotionSolver
{
	// Private data

		//- Reference point field
		pointField points0_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		displacementFvMotionSolver
		(
			const displacementFvMotionSolver&
		);

		//- Disallow default bitwise assignment
		void operator=(const displacementFvMotionSolver&);


public:

	//- Runtime type information
	TypeName("displacementInterpolation");


	// Constructors

		//- Construct from polyMesh and data stream
		displacementFvMotionSolver
		(
			const polyMesh&,
			Istream& msDataUnused
		);


	// Destructor

		virtual ~displacementFvMotionSolver();


	// Member Functions

		//- Return reference to the reference field
		const pointField& points0() const
		{
			return points0_;
		}

		//- Update topology
		virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
