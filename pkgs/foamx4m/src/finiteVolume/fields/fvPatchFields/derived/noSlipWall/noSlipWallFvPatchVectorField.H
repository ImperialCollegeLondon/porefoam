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
	Foam::noSlipWallFvPatchVectorField

Description
	Foam::noSlipWallFvPatchVectorField

SourceFiles
	noSlipWallFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef noSlipWallFvPatchVectorField_H
#define noSlipWallFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class noSlipWallFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{

public:

	//- Runtime type information
	TypeName("noSlipWall");


	// Constructors

		//- Construct from patch and internal field
		noSlipWallFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		noSlipWallFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping the given noSlipWallFvPatchVectorField
		//  onto a new patch
		noSlipWallFvPatchVectorField
		(
			const noSlipWallFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		noSlipWallFvPatchVectorField
		(
			const noSlipWallFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new noSlipWallFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		noSlipWallFvPatchVectorField
		(
			const noSlipWallFvPatchVectorField&,
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
				new noSlipWallFvPatchVectorField(*this, iF)
			);
		}


	// Member functions

			//- Return the matrix diagonal coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			virtual tmp<vectorField> gradientInternalCoeffs() const;

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			virtual tmp<vectorField> gradientBoundaryCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
