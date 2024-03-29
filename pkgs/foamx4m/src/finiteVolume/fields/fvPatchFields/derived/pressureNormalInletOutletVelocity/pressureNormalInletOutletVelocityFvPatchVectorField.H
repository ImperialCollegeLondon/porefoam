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
	Foam::pressureNormalInletOutletVelocityFvPatchVectorField

Description
	Velocity inlet/outlet boundary condition for patches where the pressure is
	specified.  zero-gradient is applied for outflow (as defined by the flux)
	and for inflow the velocity is obtained from the flux with a direction
	normal to the patch faces.

SourceFiles
	pressureNormalInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureNormalInletOutletVelocityFvPatchVectorField_H
#define pressureNormalInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class pressureNormalInletOutletVelocityFvPatchVectorField
:
	public mixedFvPatchVectorField
{
	// Private data

		//- Name of the flux field
		word phiName_;

		//- Name of the density field
		word rhoName_;


public:

	//- Runtime type information
	TypeName("pressureNormalInletOutletVelocity");


	// Constructors

		//- Construct from patch and internal field
		pressureNormalInletOutletVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		pressureNormalInletOutletVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  pressureNormalInletOutletVelocityFvPatchVectorField
		//  onto a new patch
		pressureNormalInletOutletVelocityFvPatchVectorField
		(
			const pressureNormalInletOutletVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		pressureNormalInletOutletVelocityFvPatchVectorField
		(
			const pressureNormalInletOutletVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new pressureNormalInletOutletVelocityFvPatchVectorField
				(
					*this
				)
			);
		}

		//- Construct as copy setting internal field reference
		pressureNormalInletOutletVelocityFvPatchVectorField
		(
			const pressureNormalInletOutletVelocityFvPatchVectorField&,
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
				new pressureNormalInletOutletVelocityFvPatchVectorField
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
