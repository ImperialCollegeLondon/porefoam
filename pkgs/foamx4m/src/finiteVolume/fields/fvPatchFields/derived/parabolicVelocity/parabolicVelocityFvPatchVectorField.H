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
	parabolicVelocityFvPatchVectorField

Description
	Boundary condition specifies a parabolic velocity inlet profile
	(fixed value), given maximum velocity value (peak of the parabola),
	flow direction n and direction of the parabolic coordinate y

SourceFiles
	parabolicVelocityFvPatchVectorField.C

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef parabolicVelocityFvPatchVectorField_H
#define parabolicVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class parabolicVelocityFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	// Private data

		//- Peak velocity magnitude
		scalar maxValue_;

		//- Flow direction
		vector n_;

		//- Direction of the y-coordinate
		vector y_;


public:

	//- Runtime type information
	TypeName("parabolicVelocity");


	// Constructors

		//- Construct from patch and internal field
		parabolicVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		parabolicVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given parabolicVelocityFvPatchVectorField
		//  onto a new patch
		parabolicVelocityFvPatchVectorField
		(
			const parabolicVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new parabolicVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		parabolicVelocityFvPatchVectorField
		(
			const parabolicVelocityFvPatchVectorField&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchVectorField> clone
		(
			const DimensionedField<vector, volMesh>& iF
		) const
		{
			return tmp<fvPatchVectorField>
			(
				new parabolicVelocityFvPatchVectorField(*this, iF)
			);
		}


	// Member functions

		//- Return max value
		scalar& maxValue()
		{
			return maxValue_;
		}

		//- Return flow direction
		vector& n()
		{
			return n_;
		}

		//- Return y direction
		vector& y()
		{
			return y_;
		}

		//- Update coefficients
		virtual void updateCoeffs();

		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
