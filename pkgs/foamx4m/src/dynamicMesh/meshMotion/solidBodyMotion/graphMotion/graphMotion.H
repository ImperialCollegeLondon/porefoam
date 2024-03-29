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
	Foam::solidBodyMotionFunction::graphMotion

Description
	Prescribed translational and rotational motion given graphs for surge,
	sway, heave, roll, pitch and yaw. Hence, the motion (parts of trajectory) is
	given in graphs, not the velocity.

SourceFiles
	graphMotion.C

Author
	Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef graphMotion_H
#define graphMotion_H

#include "solidBodyMotionFunction.H"
#include "graph.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{


class graphMotion
:
	public solidBodyMotionFunction
{
	// Private data

		//- Centre of gravity
		point origin_;


		// Prescribed translation members

			// - Graph for surge translational motion (x-translation)
			graph surge_;

			// - Graph for sway translational motion (y-translation)
			graph sway_;

			// - Graph for heave translational motion (z-translation)
			graph heave_;


		// Prescribed rotation members

			// - Graph for roll rotational motion (x-rotation)
			graph roll_;

			// - Graph for pitch rotational motion (y-rotation)
			graph pitch_;

			// - Graph for yaw rotational motion (z-rotation)
			graph yaw_;

			//- Is the rotational velocity given in degrees/sec
			Switch inDegrees_;


		//- Member functions

			//- Calculate transformation septernion
			septernion calcTransformation(const scalar t) const;


public:

	//- Runtime type information
	TypeName("graphMotion");


	// Constructors

		//- Construct from components
		graphMotion
		(
			const dictionary& SBMFCoeffs,
			const Time& runTime
		);


	//- Destructor
	virtual ~graphMotion()
	{}


	// Member Functions

		//- Return the solid-body motion transformation septernion
		virtual septernion transformation() const;

		//- Return the solid-body motion velocity
		virtual septernion velocity() const;

		//- Update properties from given dictionary
		virtual bool read(const dictionary& SBMFCoeffs);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
