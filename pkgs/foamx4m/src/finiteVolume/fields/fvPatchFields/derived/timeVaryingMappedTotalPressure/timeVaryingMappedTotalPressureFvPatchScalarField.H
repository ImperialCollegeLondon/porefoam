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
	Foam::timeVaryingMappedTotalPressureFvPatchScalarField

Description
	Foam::timeVaryingMappedTotalPressureFvPatchScalarField

SourceFiles
	timeVaryingMappedTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedTotalPressureFvPatchScalarField_H
#define timeVaryingMappedTotalPressureFvPatchScalarField_H

#include "timeVaryingMappedFixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class timeVaryingMappedTotalPressureFvPatchScalarField
:
	public timeVaryingMappedFixedValueFvPatchScalarField
{
	// Private data

		//- Name of the velocity field
		word UName_;

		//- Name of the flux transporting the field
		word phiName_;

		//- Name of the density field used to normalise the mass flux
		//  if neccessary
		word rhoName_;

		//- Rho value used to normalise the total pressure for incompressible flow
		scalar rho_;

		//- Name of the compressibility field used to calculate the wave speed
		word psiName_;

		//- Heat capacity ratio
		scalar gamma_;


public:

	//- Runtime type information
	TypeName("timeVaryingMappedTotalPressure");


	// Constructors

		//- Construct from patch and internal field
		timeVaryingMappedTotalPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		timeVaryingMappedTotalPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given timeVaryingMappedTotalPressureFvPatchScalarField
		//  onto a new patch
		timeVaryingMappedTotalPressureFvPatchScalarField
		(
			const timeVaryingMappedTotalPressureFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		timeVaryingMappedTotalPressureFvPatchScalarField
		(
			const timeVaryingMappedTotalPressureFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				new timeVaryingMappedTotalPressureFvPatchScalarField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		timeVaryingMappedTotalPressureFvPatchScalarField
		(
			const timeVaryingMappedTotalPressureFvPatchScalarField&,
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
				new timeVaryingMappedTotalPressureFvPatchScalarField(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return the name of the velocity field
			const word& UName() const
			{
				return UName_;
			}

			//- Return reference to the name of the velocity field
			//  to allow adjustment
			word& UName()
			{
				return UName_;
			}

			//- Return the heat capacity ratio
			scalar gamma() const
			{
				return gamma_;
			}

			//- Return reference to the heat capacity ratio to allow adjustment
			scalar& gamma()
			{
				return gamma_;
			}


		// Mapping functions

			//- Map (and resize as needed) from self given a mapping object
			virtual void autoMap
			(
				const fvPatchFieldMapper&
			);

			//- Reverse map the given fvPatchField onto this fvPatchField
			virtual void rmap
			(
				const fvPatchScalarField&,
				const labelList&
			);


		// Evaluation functions

			//- Update the coefficients associated with the patch field
			//  using the given patch velocity field
			virtual void updateCoeffs(const vectorField& Up);

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
