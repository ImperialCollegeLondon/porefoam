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
	Foam::fluxCorrectedVelocityFvPatchVectorField

Description
	Velocity outlet boundary condition for patches where the pressure is
	specified. The outflow velocity is obtained by "zeroGradient" and then
	corrected from the flux.  If reverse flow is possible or expected use
	the "pressureInletOutletVelocityFvPatchVectorField" BC instead.

SourceFiles
	fluxCorrectedVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fluxCorrectedVelocityFvPatchVectorField_H
#define fluxCorrectedVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class fluxCorrectedVelocityFvPatchVectorField
:
	public zeroGradientFvPatchVectorField
{
	// Private data

		//- Name of flux field
		word phiName_;

		//- Name of density field
		word rhoName_;


public:

	//- Runtime type information
	TypeName("fluxCorrectedVelocity");


	// Constructors

		//- Construct from patch and internal field
		fluxCorrectedVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		fluxCorrectedVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given fluxCorrectedVelocityFvPatchVectorField
		//  onto a new patch
		fluxCorrectedVelocityFvPatchVectorField
		(
			const fluxCorrectedVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new fluxCorrectedVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		fluxCorrectedVelocityFvPatchVectorField
		(
			const fluxCorrectedVelocityFvPatchVectorField&,
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
				new fluxCorrectedVelocityFvPatchVectorField(*this, iF)
			);
		}


	// Member functions

		// Evaluation functions

			//- Evaluate the patch field
			virtual void evaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			);


		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
