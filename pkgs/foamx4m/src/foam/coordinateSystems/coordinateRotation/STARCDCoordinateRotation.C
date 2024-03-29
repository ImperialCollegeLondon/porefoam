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

#include "STARCDCoordinateRotation.H"

#include "Switch.H"
#include "mathematicalConstants.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(STARCDCoordinateRotation, 0);
	addToRunTimeSelectionTable
	(
		coordinateRotation,
		STARCDCoordinateRotation,
		dictionary
	);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::STARCDCoordinateRotation::calcTransform
(
	const scalar rotZ,
	const scalar rotX,
	const scalar rotY,
	const bool inDegrees
)
{
	scalar x = rotX;
	scalar y = rotY;
	scalar z = rotZ;

	if (inDegrees)
	{
		x *= mathematicalConstant::pi/180.0;
		y *= mathematicalConstant::pi/180.0;
		z *= mathematicalConstant::pi/180.0;
	}

	tensor::operator=
	(
		tensor
		(
			cos(y)*cos(z) - sin(x)*sin(y)*sin(z),
			-cos(x)*sin(z),
			sin(x)*cos(y)*sin(z) + sin(y)*cos(z),

			cos(y)*sin(z) + sin(x)*sin(y)*cos(z),
			cos(x)*cos(z),
			sin(y)*sin(z) - sin(x)*cos(y)*cos(z),

			-cos(x)*sin(y),
			sin(x),
			cos(x)*cos(y)
		)
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::STARCDCoordinateRotation::STARCDCoordinateRotation()
:
	coordinateRotation()
{}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
	const vector& rotZrotXrotY,
	const bool inDegrees
)
:
	coordinateRotation()
{
	calcTransform
	(
		rotZrotXrotY.component(vector::X),
		rotZrotXrotY.component(vector::Y),
		rotZrotXrotY.component(vector::Z),
		inDegrees
	);
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
	const scalar rotZ,
	const scalar rotX,
	const scalar rotY,
	const bool inDegrees
)
:
	coordinateRotation()
{
	calcTransform(rotZ, rotX, rotY, inDegrees);
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
	const dictionary& dict
)
:
	coordinateRotation()
{
	vector rotation(dict.lookup("rotation"));

	calcTransform
	(
		rotation.component(vector::X),
		rotation.component(vector::Y),
		rotation.component(vector::Z),
		dict.lookupOrDefault<Switch>("degrees", true)
	);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
