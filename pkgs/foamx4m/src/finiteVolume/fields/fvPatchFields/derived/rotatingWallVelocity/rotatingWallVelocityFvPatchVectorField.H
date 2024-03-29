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
	Foam::rotatingWallVelocityFvPatchVectorField

Description
	Foam::rotatingWallVelocityFvPatchVectorField

SourceFiles
	rotatingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingWallVelocityFvPatchVectorField_H
#define rotatingWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class rotatingWallVelocityFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	// Private data

		//- Origin of the rotation
		vector origin_;

		//- Axis of the rotation
		vector axis_;

		//- Rotational speed
		scalar omega_;


public:

	//- Runtime type information
	TypeName("rotatingWallVelocity");


	// Constructors

		//- Construct from patch and internal field
		rotatingWallVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		rotatingWallVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given rotatingWallVelocityFvPatchVectorField
		//  onto a new patch
		rotatingWallVelocityFvPatchVectorField
		(
			const rotatingWallVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		rotatingWallVelocityFvPatchVectorField
		(
			const rotatingWallVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new rotatingWallVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		rotatingWallVelocityFvPatchVectorField
		(
			const rotatingWallVelocityFvPatchVectorField&,
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
				new rotatingWallVelocityFvPatchVectorField(*this, iF)
			);
		}



	// Member functions

		// Access functions

			//- Return the origin of the rotation
			const vector& origin() const
			{
				return origin_;
			}

			//- Return the axis of the rotation
			const vector& axis() const
			{
				return axis_;
			}

			//- Return the rotational speed
			const scalar& omega() const
			{
				return omega_;
			}

			//- Return non-const access to the origin of the rotation
			vector& origin()
			{
				return origin_;
			}

			//- Return non-const access to the axis of the rotation
			vector& axis()
			{
				return axis_;
			}

			//- Return non-const access to the rotational speed
			scalar& omega()
			{
				return omega_;
			}


		//- Update the coefficients associated with the patch field
		virtual void updateCoeffs();

		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
