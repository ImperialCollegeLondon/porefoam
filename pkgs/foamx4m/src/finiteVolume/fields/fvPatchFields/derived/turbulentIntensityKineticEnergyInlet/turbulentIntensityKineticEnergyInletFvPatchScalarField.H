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
	Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField

Description
	Calculate turbulent kinetic energy from the intensity provided as a
	fraction of the mean velocity

	Example of the boundary condition specification:
	\verbatim
		inlet
		{
			type		turbulentIntensityKineticEnergyInlet;
			U		   U;			  // name of velocity field
			intensity   0.05;		   // 5% turbulence
			value	   uniform 1;	  // placeholder
		}
	\endverbatim

SourceFiles
	turbulentIntensityKineticEnergyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentIntensityKineticEnergyInletFvPatchScalarField_H
#define turbulentIntensityKineticEnergyInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class turbulentIntensityKineticEnergyInletFvPatchScalarField
:
	public inletOutletFvPatchScalarField
{
	// Private data

		//- Turbulent intensity as fraction of mean velocity
		scalar intensity_;

		//- Name of the velocity field
		word UName_;

		//- Name of the flux field
		word phiName_;


public:

	//- Runtime type information
	TypeName("turbulentIntensityKineticEnergyInlet");


	// Constructors

		//- Construct from patch and internal field
		turbulentIntensityKineticEnergyInletFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		turbulentIntensityKineticEnergyInletFvPatchScalarField
		(
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  turbulentIntensityKineticEnergyInletFvPatchScalarField
		//  onto a new patch
		turbulentIntensityKineticEnergyInletFvPatchScalarField
		(
			const turbulentIntensityKineticEnergyInletFvPatchScalarField&,
			const fvPatch&,
			const DimensionedField<scalar, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		turbulentIntensityKineticEnergyInletFvPatchScalarField
		(
			const turbulentIntensityKineticEnergyInletFvPatchScalarField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchScalarField> clone() const
		{
			return tmp<fvPatchScalarField>
			(
				new turbulentIntensityKineticEnergyInletFvPatchScalarField
				(
					*this
				)
			);
		}

		//- Construct as copy setting internal field reference
		turbulentIntensityKineticEnergyInletFvPatchScalarField
		(
			const turbulentIntensityKineticEnergyInletFvPatchScalarField&,
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
				new turbulentIntensityKineticEnergyInletFvPatchScalarField
				(
					*this,
					iF
				)
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
