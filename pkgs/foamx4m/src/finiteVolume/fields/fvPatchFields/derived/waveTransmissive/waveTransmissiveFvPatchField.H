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
	Foam::waveTransmissiveFvPatchField

Description
	Foam::waveTransmissiveFvPatchField

SourceFiles
	waveTransmissiveFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef waveTransmissiveFvPatchField_H
#define waveTransmissiveFvPatchField_H

#include "advectiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class waveTransmissiveFvPatchField
:
	public advectiveFvPatchField<Type>
{

	// Private data

		//- Name of the compressibility field used to calculate the wave speed
		word psiName_;

		//- Name of velocity field used to calculate supercritical condition
		word UName_;

		//- Heat capacity ratio
		scalar gamma_;


public:

	//- Runtime type information
	TypeName("waveTransmissive");


	// Constructors

		//- Construct from patch and internal field
		waveTransmissiveFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		waveTransmissiveFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given waveTransmissiveFvPatchField
		//  onto a new patch
		waveTransmissiveFvPatchField
		(
			const waveTransmissiveFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		waveTransmissiveFvPatchField
		(
			const waveTransmissiveFvPatchField&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new waveTransmissiveFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		waveTransmissiveFvPatchField
		(
			const waveTransmissiveFvPatchField&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchField<Type> > clone
		(
			const DimensionedField<Type, volMesh>& iF
		) const
		{
			return tmp<fvPatchField<Type> >
			(
				new waveTransmissiveFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return the heat capacity ratio
			scalar gamma() const
			{
				return gamma_;
			}

			//- Return reference to the heat capacity ratio to allow adjustment
			scalar& gamma()
			{
				return gamma_;
			}


		// Evaluation functions

			//- Calculate and return the advection speed at the boundary
			virtual tmp<scalarField> advectionSpeed() const;

			//- Calculate and return the supercritical switch at the boundary
			//  Supercritical = 1 converts the outlet to zeroGradient
			//  Supercritical = 0 no correction
			virtual tmp<scalarField> supercritical() const;


		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "waveTransmissiveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
