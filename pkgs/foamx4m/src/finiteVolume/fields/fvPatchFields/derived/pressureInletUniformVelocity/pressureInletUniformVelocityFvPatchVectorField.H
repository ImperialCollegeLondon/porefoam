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
	Foam::pressureInletUniformVelocityFvPatchVectorField

Description
	Velocity inlet boundary condition for patches where the pressure is
	specified. The uniform inflow velocity is obtained by averaging the flux
	over the patch and apply it in the direction normal to the patch faces.

SourceFiles
	pressureInletUniformVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInletUniformVelocityFvPatchVectorField_H
#define pressureInletUniformVelocityFvPatchVectorField_H

#include "pressureInletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class pressureInletUniformVelocityFvPatchVectorField
:
	public pressureInletVelocityFvPatchVectorField
{

public:

	//- Runtime type information
	TypeName("pressureInletUniformVelocity");


	// Constructors

		//- Construct from patch and internal field
		pressureInletUniformVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		pressureInletUniformVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given
		//  pressureInletUniformVelocityFvPatchVectorField onto a new patch
		pressureInletUniformVelocityFvPatchVectorField
		(
			const pressureInletUniformVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		pressureInletUniformVelocityFvPatchVectorField
		(
			const pressureInletUniformVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new pressureInletUniformVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		pressureInletUniformVelocityFvPatchVectorField
		(
			const pressureInletUniformVelocityFvPatchVectorField&,
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
				new pressureInletUniformVelocityFvPatchVectorField(*this, iF)
			);
		}


	// Member functions

		//- Update the coefficients associated with the patch field
		virtual void updateCoeffs();


	// Member operators

		virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
