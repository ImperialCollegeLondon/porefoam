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
	Foam::fixedFluxPressureFvPatchScalarField

Description
	Foam::fixedFluxPressureFvPatchScalarField

SourceFiles
	fixedFluxPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedFluxPressureFvPatchScalarFields_H
#define fixedFluxPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class fixedFluxPressureFvPatchScalarField
:
	public fixedGradientFvPatchScalarField
{
	// Private data

		//- Name of the velocity field
		word UName_;

		//- Name of the flux transporting the field
		word phiName_;

		//- Name of the density field used to normalise the mass flux
		//  if neccessary
		word rhoName_;

		//- Is the pressure adjoint, i.e. has the opposite sign
		word rAUName_;

		//- Is the pressure adjoint, i.e. has the opposite sign
		Switch adjoint_;


public:

	//- Runtime type information
	TypeName("fixedFluxPressure");


	// Constructors

		//- Construct from patch and internal field
		fixedFluxPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		fixedFluxPressureFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given fixedFluxPressureFvPatchScalarField onto
		//  a new patch
		fixedFluxPressureFvPatchScalarField
		(
			const fixedFluxPressureFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		fixedFluxPressureFvPatchScalarField
		(
			const fixedFluxPressureFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				new fixedFluxPressureFvPatchScalarField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		fixedFluxPressureFvPatchScalarField
		(
			const fixedFluxPressureFvPatchScalarField&,
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
				new fixedFluxPressureFvPatchScalarField(*this, iF)
			);
		}


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
