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
	Foam::STARCDCoordinateRotation

Description
	A coordinateRotation defined by the STAR-CD convention.

	The 3 rotations are defined in the STAR-CD convention
	(around Z, around X' and around Y'').
	The order of the parameter arguments matches this rotation order.

	- the rotation angles are in degrees, unless otherwise explictly specified:

	@verbatim
		coordinateRotation
		{
			type        STARCDRotation;
			degrees     false;
			rotation    (0 0 3.141592654);
		}
	@endverbatim

\*---------------------------------------------------------------------------*/

#ifndef STARCDCoordinateRotation_H
#define STARCDCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class STARCDCoordinateRotation
:
	public coordinateRotation
{
	// Private member functions

		//- Calculate transformation tensor
		void calcTransform
		(
			const scalar rotZ,
			const scalar rotX,
			const scalar rotY,
			const bool inDegrees = true
		);


public:

	//- Runtime type information
	TypeName("STARCDRotation");


	// Constructors

		//- Construct null
		STARCDCoordinateRotation();

		//- Construct from rotation vector
		STARCDCoordinateRotation
		(
			const vector& rotZrotXrotY,
			const bool inDegrees = true
		);

		//- Construct from components of rotation vector
		STARCDCoordinateRotation
		(
			const scalar rotZ,
			const scalar rotX,
			const scalar rotY,
			const bool inDegrees = true
		);

		//- Construct from dictionary
		STARCDCoordinateRotation(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
