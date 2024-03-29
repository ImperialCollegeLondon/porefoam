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
	Foam::freestreamPressureFvPatchScalarField

Description
	Foam::freestreamPressureFvPatchScalarField

SourceFiles
	freestreamPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef freestreamPressureFvPatchScalarFields_H
#define freestreamPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class freestreamPressureFvPatchScalarField
:
	public zeroGradientFvPatchScalarField
{

public:

	//- Runtime type information
	TypeName("freestreamPressure");


	// Constructors

		//- Construct from patch and internal field
		freestreamPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		freestreamPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given freestreamPressureFvPatchScalarField onto
		//  a new patch
		freestreamPressureFvPatchScalarField
		(
			const freestreamPressureFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		freestreamPressureFvPatchScalarField
		(
			const freestreamPressureFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				new freestreamPressureFvPatchScalarField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		freestreamPressureFvPatchScalarField
		(
			const freestreamPressureFvPatchScalarField&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchScalarField> clone
		(
			const DimensionedField<scalar, volMesh>& iF
		) const
		{
			return tmp<fvPatchScalarField>
			(
				new freestreamPressureFvPatchScalarField(*this, iF)
			);
		}


	// Member functions

		// Evaluation functions

			//- Update the coefficients associated with the patch field
			virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
