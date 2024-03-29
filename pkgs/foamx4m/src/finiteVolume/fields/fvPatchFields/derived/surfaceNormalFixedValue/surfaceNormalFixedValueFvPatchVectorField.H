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
	Foam::surfaceNormalFixedValueFvPatchVectorField

Description
	Describes a surface normal vector boundary condition by its magnitude.
	Note: The value is positive for outward-pointing vectors

SourceFiles
	surfaceNormalFixedValueFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceNormalFixedValueFvPatchVectorField_H
#define surfaceNormalFixedValueFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class surfaceNormalFixedValueFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	// Private data

		//- Surface-normal velocity value
		scalarField refValue_;


public:

	//- Runtime type information
	TypeName("surfaceNormalFixedValue");


	// Constructors

		//- Construct from patch and internal field
		surfaceNormalFixedValueFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		surfaceNormalFixedValueFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  surfaceNormalFixedValueFvPatchVectorField
		//  onto a new patch
		surfaceNormalFixedValueFvPatchVectorField
		(
			const surfaceNormalFixedValueFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		surfaceNormalFixedValueFvPatchVectorField
		(
			const surfaceNormalFixedValueFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new surfaceNormalFixedValueFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		surfaceNormalFixedValueFvPatchVectorField
		(
			const surfaceNormalFixedValueFvPatchVectorField&,
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
				new surfaceNormalFixedValueFvPatchVectorField
				(
					*this,
					iF
				)
			);
		}



	// Member functions

		// Mapping functions

			//- Map (and resize as needed) from self given a mapping object
			virtual void autoMap
			(
				const fvPatchFieldMapper&
			);

			//- Reverse map the given fvPatchField onto this fvPatchField
			virtual void rmap
			(
				const fvPatchVectorField&,
				const labelList&
			);


	// Member functions

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
