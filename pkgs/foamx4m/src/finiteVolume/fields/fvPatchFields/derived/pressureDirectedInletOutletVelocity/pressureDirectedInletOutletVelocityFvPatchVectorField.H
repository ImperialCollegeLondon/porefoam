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
	Foam::pressureDirectedInletOutletVelocityFvPatchVectorField

Description
	Velocity inlet/outlet boundary condition for pressure boundary where the
	pressure is specified.  zero-gradient is applied for outflow (as defined
	by the flux) and for inflow the velocity is obtained from the flux with
	the specified `inletDirection'.

SourceFiles
	pressureDirectedInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureDirectedInletOutletVelocityFvPatchVectorField_H
#define pressureDirectedInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class pressureDirectedInletOutletVelocityFvPatchVectorField
:
	public mixedFvPatchVectorField
{
	// Private data

		//- Name of flux field
		word phiName_;

		//- Name of density field
		word rhoName_;

		//- Inlet direction
		vectorField inletDir_;


public:

	//- Runtime type information
	TypeName("pressureDirectedInletOutletVelocity");


	// Constructors

		//- Construct from patch and internal field
		pressureDirectedInletOutletVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		pressureDirectedInletOutletVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  pressureDirectedInletOutletVelocityFvPatchVectorField
		//  onto a new patch
		pressureDirectedInletOutletVelocityFvPatchVectorField
		(
			const pressureDirectedInletOutletVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		pressureDirectedInletOutletVelocityFvPatchVectorField
		(
			const pressureDirectedInletOutletVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new pressureDirectedInletOutletVelocityFvPatchVectorField
				(
					*this
				)
			);
		}

		//- Construct as copy setting internal field reference
		pressureDirectedInletOutletVelocityFvPatchVectorField
		(
			const pressureDirectedInletOutletVelocityFvPatchVectorField&,
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
				new pressureDirectedInletOutletVelocityFvPatchVectorField
				(
					*this,
					iF
				)
			);
		}


	// Member functions

		// Access

			//- Return the name of rho
			const word& rhoName() const
			{
				return rhoName_;
			}

			//- Return reference to the name of rho to allow adjustment
			word& rhoName()
			{
				return rhoName_;
			}

			//- Return the name of phi
			const word& phiName() const
			{
				return phiName_;
			}

			//- Return reference to the name of phi to allow adjustment
			word& phiName()
			{
				return phiName_;
			}

			//- Return inlet direction
			const vectorField& inletDir() const
			{
				return inletDir_;
			}

			//- Return access to inlet direction
			vectorField& inletDir()
			{
				return inletDir_;
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
				const fvPatchVectorField&,
				const labelList&
			);


		//- Update the coefficients associated with the patch field
		virtual void updateCoeffs();

		//- Write
		virtual void write(Ostream&) const;


	// Member operators

		virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
