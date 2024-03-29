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
	Foam::uniformDensityHydrostaticPressureFvPatchScalarField

Description
	Hydrostatic pressure boundary condition calculated as

	pRefValue + rho*g.(x - pRefPoint)

	where rho is provided and assumed uniform.

SourceFiles
	uniformDensityHydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDensityHydrostaticPressureFvPatchScalarField_H
#define uniformDensityHydrostaticPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class uniformDensityHydrostaticPressureFvPatchScalarField
:
	public fixedValueFvPatchScalarField
{
	// Private data

		//- Constant density in the far-field
		scalar rho_;

		//- Reference pressure
		scalar pRefValue_;

		//- Reference pressure location
		vector pRefPoint_;


public:

	//- Runtime type information
	TypeName("uniformDensityHydrostaticPressure");


	// Constructors

		//- Construct from patch and internal field
		uniformDensityHydrostaticPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		uniformDensityHydrostaticPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  uniformDensityHydrostaticPressureFvPatchScalarField onto a new patch
		uniformDensityHydrostaticPressureFvPatchScalarField
		(
			const uniformDensityHydrostaticPressureFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		uniformDensityHydrostaticPressureFvPatchScalarField
		(
			const uniformDensityHydrostaticPressureFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				new uniformDensityHydrostaticPressureFvPatchScalarField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		uniformDensityHydrostaticPressureFvPatchScalarField
		(
			const uniformDensityHydrostaticPressureFvPatchScalarField&,
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
				new uniformDensityHydrostaticPressureFvPatchScalarField
				(
					*this,
					iF
				)
			);
		}


	// Member functions

		// Access

			//- Return the constant density in the far-field
			scalar rho() const
			{
				return rho_;
			}

			//- Return reference to the constant density in the far-field
			//  to allow adjustment
			scalar& rho()
			{
				return rho_;
			}

			//- Return the reference pressure
			scalar pRefValue() const
			{
				return pRefValue_;
			}

			//- Return reference to the reference pressure to allow adjustment
			scalar& pRefValue()
			{
				return pRefValue_;
			}

			//- Return the pressure reference location
			const vector& pRefPoint() const
			{
				return pRefPoint_;
			}

			//- Return reference to the pressure reference location
			//  to allow adjustment
			vector& pRefPoint()
			{
				return pRefPoint_;
			}


		// Evaluation functions

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
