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
	Foam::uniformDiffusivity

Description
	Uniform uniform finite volume mesh motion diffusivity.

SourceFiles
	uniformDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDiffusivity_H
#define uniformDiffusivity_H

#include "motionDiffusivity.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class uniformDiffusivity
:
	public motionDiffusivity
{

protected:

	// Protected data

		surfaceScalarField faceDiffusivity_;


private:

	// Private Member Functions

		//- Disallow default bitwise copy construct
		uniformDiffusivity(const uniformDiffusivity&);

		//- Disallow default bitwise assignment
		void operator=(const uniformDiffusivity&);


public:

	//- Runtime type information
	TypeName("uniform");


	// Constructors

		//- Construct for the given fvMotionSolver and data Istream
		uniformDiffusivity
		(
			const fvMotionSolver& mSolver,
			Istream& mdData
		);


	// Destructor

		virtual ~uniformDiffusivity();


	// Member Functions

		//- Return diffusivity field
		virtual tmp<surfaceScalarField> operator()() const
		{
			return faceDiffusivity_;
		}

		//- Do not correct the motion diffusivity
		virtual void correct()
		{}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
