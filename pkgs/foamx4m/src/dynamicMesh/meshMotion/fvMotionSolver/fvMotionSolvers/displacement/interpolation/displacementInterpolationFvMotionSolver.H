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
	Foam::displacementInterpolationFvMotionSolver

Description
	Mesh motion solver for an fvMesh.

	Scales inbetween motion prescribed on faceZones. Works out per point
	the distance between the bounding face zones (in all three directions)
	at the start and then every time step
	- moves the faceZones based on tables
	- interpolates the displacement of all points based on the
	  faceZone motion.

	Tables are in the @a constant/tables directory.

Note
	could be a motionSolver - does not use any fvMesh structure.

SourceFiles
	displacementInterpolationFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementInterpolationFvMotionSolver_H
#define displacementInterpolationFvMotionSolver_H

#include "displacementFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class displacementInterpolationFvMotionSolver
:
	public displacementFvMotionSolver
{
	// Private data

		//- Additional settings for motion solver
		dictionary dynamicMeshCoeffs_;


		// Face zone information (note: could pack these to only contain
		// used zones)

			//- Interpolation table. From faceZone to times.
			List<scalarField> times_;

			//- Interpolation table. From faceZone to displacements.
			List<vectorField> displacements_;

		// Range information.

			//- Per direction, per range the index of the lower
			//  faceZone
			FixedList<labelList, 3> rangeToZone_;

			//- Per direction, per range the points that are in it
			FixedList<labelListList, 3> rangeToPoints_;

			//- Per direction, per range the weight of the points relative
			//  to this and the next range.
			FixedList<List<scalarField>, 3> rangeToWeights_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		displacementInterpolationFvMotionSolver
		(
			const displacementInterpolationFvMotionSolver&
		);

		//- Disallow default bitwise assignment
		void operator=(const displacementInterpolationFvMotionSolver&);


public:

	//- Runtime type information
	TypeName("displacementInterpolation");


	// Constructors

		//- Construct from polyMesh and data stream
		displacementInterpolationFvMotionSolver
		(
			const polyMesh&,
			Istream& msDataUnused
		);


	// Destructor

		~displacementInterpolationFvMotionSolver();


	// Member Functions

		//- Return point location obtained from the current motion field
		virtual tmp<pointField> curPoints() const;

		//- Solve for motion
		virtual void solve()
		{}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
