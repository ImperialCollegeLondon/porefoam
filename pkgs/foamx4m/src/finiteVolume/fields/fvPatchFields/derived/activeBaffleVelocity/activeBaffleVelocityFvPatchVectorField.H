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
	Foam::activeBaffleVelocityFvPatchVectorField

Description
	Boundary condition that modifies mesh areas based on velocity.

SourceFiles
	activeBaffleVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef activeBaffleVelocityFvPatchVectorField_H
#define activeBaffleVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class activeBaffleVelocityFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	// Private data

		//- Name of the pressure field used to calculate the force
		//  on the active baffle
		word pName_;

		//- Name of the cyclic patch used when the active baffle is open
		word cyclicPatchName_;

		//- Index of the cyclic patch used when the active baffle is open
		label cyclicPatchLabel_;

		//- Orientation (1 or -1) of the active baffle patch.
		//  Used to change the direction of opening without the need for
		//  reordering the patch faces
		label orientation_;

		//- Initial wall patch areas
		vectorField initWallSf_;

		//- Initial cyclic patch areas
		vectorField initCyclicSf_;

		//- Current fraction of the active baffle which is open
		scalar openFraction_;

		//- Time taken for the active baffle to open
		scalar openingTime_;

		//- Maximum fractional change to the active baffle openness
		//  per time-step
		scalar maxOpenFractionDelta_;

		label curTimeIndex_;


public:

	//- Runtime type information
	TypeName("activeBaffleVelocity");


	// Constructors

		//- Construct from patch and internal field
		activeBaffleVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		activeBaffleVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given activeBaffleVelocityFvPatchVectorField
		//  onto a new patch
		activeBaffleVelocityFvPatchVectorField
		(
			const activeBaffleVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		activeBaffleVelocityFvPatchVectorField
		(
			const activeBaffleVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new activeBaffleVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		activeBaffleVelocityFvPatchVectorField
		(
			const activeBaffleVelocityFvPatchVectorField&,
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
				new activeBaffleVelocityFvPatchVectorField(*this, iF)
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
