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
	Foam::SRFVelocityFvPatchVectorField

Description
	Velocity patch to be used with SRF model

SourceFiles
	SRFVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef SRFVelocityFvPatchVectorField_H
#define SRFVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class SRFVelocityFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	// Private data

		//- Is the supplied inlet value relative to the SRF
		Switch relative_;

		//- Inlet value
		vectorField inletValue_;


public:

	//- Runtime type information
	TypeName("SRFVelocity");


	// Constructors

		//- Construct from patch and internal field
		SRFVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		SRFVelocityFvPatchVectorField
		(
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given SRFVelocityFvPatchVectorField
		//  onto a new patch
		SRFVelocityFvPatchVectorField
		(
			const SRFVelocityFvPatchVectorField&,
			const fvPatch&,
			const DimensionedField<vector, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		SRFVelocityFvPatchVectorField
		(
			const SRFVelocityFvPatchVectorField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchVectorField> clone() const
		{
			return tmp<fvPatchVectorField>
			(
				new SRFVelocityFvPatchVectorField(*this)
			);
		}

		//- Construct as copy setting internal field reference
		SRFVelocityFvPatchVectorField
		(
			const SRFVelocityFvPatchVectorField&,
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
				new SRFVelocityFvPatchVectorField(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return const access to the relative flag
			const Switch& relative() const
			{
				return relative_;
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


		// Evaluation functions

			//- Update the coefficients associated with the patch field
			virtual void updateCoeffs();


		// I-O

			//- Write
			virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
