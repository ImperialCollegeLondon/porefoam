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
	Foam::displacementSBRStressFvMotionSolver

Description
	Mesh motion solver for an fvMesh.  Based on solving the cell-centre
	solid-body rotation stress equations for the motion displacement.

SourceFiles
	displacementSBRStressFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementSBRStressFvMotionSolver_H
#define displacementSBRStressFvMotionSolver_H

#include "displacementFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;


class displacementSBRStressFvMotionSolver
:
	public displacementFvMotionSolver
{
	// Private data

		//- Point motion field
		mutable pointVectorField pointDisplacement_;

		//- Cell-centre motion field
		mutable volVectorField cellDisplacement_;

		//- Diffusivity used to control the motion
		autoPtr<motionDiffusivity> diffusivityPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		displacementSBRStressFvMotionSolver
		(
			const displacementSBRStressFvMotionSolver&
		);

		//- Disallow default bitwise assignment
		void operator=(const displacementSBRStressFvMotionSolver&);


public:

	//- Runtime type information
	TypeName("displacementSBRStress");


	// Constructors

		//- Construct from polyMesh and data stream
		displacementSBRStressFvMotionSolver
		(
			const polyMesh&,
			Istream& msDataUnused
		);


	// Destructor

		~displacementSBRStressFvMotionSolver();


	// Member Functions

		//- Return reference to the point motion displacement field
		pointVectorField& pointDisplacement()
		{
			return pointDisplacement_;
		}

		//- Return const reference to the point motion displacement field
		const pointVectorField& pointDisplacement() const
		{
			return pointDisplacement_;
		}

		//- Return reference to the cell motion displacement field
		volVectorField& cellDisplacement()
		{
			return cellDisplacement_;
		}

		//- Return const reference to the cell motion displacement field
		const volVectorField& cellDisplacement() const
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
