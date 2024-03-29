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
	Foam::totalTemperatureFvPatchScalarField

Description
	Total temperature boundary condition.  The snGrad is set to zero to stop
	back-diffusion

SourceFiles
	totalTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalTemperatureFvPatchScalarField_H
#define totalTemperatureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class totalTemperatureFvPatchScalarField
:
	public fixedValueFvPatchScalarField
{
	// Private data

		//- Name of the velocity field
		word UName_;

		//- Name of the flux transporting the field
		word phiName_;

		//- Name of the compressibility field used to calculate the wave speed
		word psiName_;

		//- Heat capacity ratio
		scalar gamma_;

		//- Total temperature
		scalarField T0_;


public:

	//- Runtime type information
	TypeName("totalTemperature");


	// Constructors

		//- Construct from patch and internal field
		totalTemperatureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		totalTemperatureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given totalTemperatureFvPatchScalarField
		//  onto a new patch
		totalTemperatureFvPatchScalarField
		(
			const totalTemperatureFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		totalTemperatureFvPatchScalarField
		(
			const totalTemperatureFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				new totalTemperatureFvPatchScalarField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		totalTemperatureFvPatchScalarField
		(
			const totalTemperatureFvPatchScalarField&,
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
				new totalTemperatureFvPatchScalarField(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return the name of the velocity field
			const word& UName() const
			{
				return UName_;
			}

			//- Return the total pressure
			const scalarField& T0() const
			{
				return T0_;
			}

			//- Return reference to the total pressure to allow adjustment
			scalarField& T0()
			{
				return T0_;
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
			virtual void updateCoeffs();

			//- Return patch-normal gradient: set to zero
			virtual tmp<scalarField> snGrad() const;

			//- Return the matrix diagonal coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			// Removed gradient contribution
			virtual tmp<scalarField> gradientInternalCoeffs() const;

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			virtual tmp<scalarField> gradientBoundaryCoeffs() const;


		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
