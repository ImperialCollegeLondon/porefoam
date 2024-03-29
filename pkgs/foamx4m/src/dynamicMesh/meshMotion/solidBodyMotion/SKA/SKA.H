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
	Foam::solidBodyMotionFunctions::SKA

Description
	Sea Keeping Analysis (SKA) 6DoF motion function.

	Obtained by interpolating tabulated data for surge (x-translation),
	sway (y-translation), heave (z-translation), roll (rotation about x),
	pitch (rotation about y) and yaw (rotation about z).

See Also
	SDA (Ship design analysis) for 3DoF motion.

SourceFiles
	SKA.C

\*---------------------------------------------------------------------------*/

#ifndef SKA_H
#define SKA_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "Vector2DTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{


class SKA
:
	public solidBodyMotionFunction
{
	// Private data

		//- Time data file name read from dictionary
		fileName timeDataFileName_;

		//- Center of gravity read from dictionary
		vector CofG_;

		//- Type used to read in the translation and rotation "vectors"
		typedef Vector2D<vector> translationRotationVectors;

		//- Field of times
		scalarField times_;

		//- Field of translation and rotation "vectors"
		Field<translationRotationVectors> values_;


	// Private Member Functions

		//- Disallow copy construct
		SKA(const SKA&);

		//- Disallow default bitwise assignment
		void operator=(const SKA&);


		//- Calculate transform given time
		septernion calcTransformation(const scalar t) const;


public:

	//- Runtime type information
	TypeName("SKA");


	// Constructors

		//- Construct from components
		SKA
		(
			const dictionary& SBMFCoeffs,
			const Time& runTime
		);


	// Destructor

		virtual ~SKA();


	// Member Functions

		//- Return the solid-body motion transformation septernion
		virtual septernion transformation() const;

		//- Return the solid-body motion velocity septernion
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
