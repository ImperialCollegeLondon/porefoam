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
	Foam::flowRateInletVelocityFvPatchVectorField

Description
	Describes a volumetric/mass flow normal vector boundary condition by its
	magnitude as an integral over its area.

	The basis of the patch (volumetric or mass) is determined by the
	dimensions of the flux, phi.
	The current density is used to correct the velocity when applying the
	mass basis.

	Example of the boundary condition specification:
	@verbatim
	inlet
	{
		type			flowRateInletVelocity;
		flowRate		0.2;		// Volumetric/mass flow rate [m3/s or kg/s]
		value		   uniform (0 0 0); // placeholder
	}
	@endverbatim

Note
	- The value is positive inwards
	- May not work correctly for transonic inlets
	- Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
	flowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletVelocityFvPatchVectorField_H
#define flowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class flowRateInletVelocityFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	// Private data

		//- Inlet integral flow rate
		scalar flowRate_;

		//- Name of the flux transporting the field
		word phiName_;

		//- Name of the density field used to normalize the mass flux
		word rhoName_;

		//- Sum areas.  Note: not correctly updated on mesh motion
		scalar gSumArea_;


public:

   //- Runtime type information
   TypeName("flowRateInletVelocity");


   // Constructors

		//- Construct from patch and internal field
		flowRateInletVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		flowRateInletVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  flowRateInletVelocityFvPatchVectorField
		//  onto a new patch
		flowRateInletVelocityFvPatchVectorField
		(
			const flowRateInletVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		flowRateInletVelocityFvPatchVectorField
		(
			const flowRateInletVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new flowRateInletVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		flowRateInletVelocityFvPatchVectorField
		(
			const flowRateInletVelocityFvPatchVectorField&,
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
				new flowRateInletVelocityFvPatchVectorField(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return the flux
			scalar flowRate() const
			{
				return flowRate_;
			}

			//- Return reference to the flux to allow adjustment
			scalar& flowRate()
			{
				return flowRate_;
			}


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
