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
	Foam::cyclicFvsPatchField

Description
	Foam::cyclicFvsPatchField

SourceFiles
	cyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFvsPatchField_H
#define cyclicFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class cyclicFvsPatchField
:
	public coupledFvsPatchField<Type>
{
	// Private data

		//- Local reference cast into the cyclic patch
		const cyclicFvPatch& cyclicPatch_;


public:

	//- Runtime type information
	TypeName(cyclicFvPatch::typeName_());


	// Constructors

		//- Construct from patch and internal field
		cyclicFvsPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, surfaceMesh>&
		);

		//- Construct from patch, internal field and dictionary
		cyclicFvsPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, surfaceMesh>&,
			const dictionary&
		);

		//- Construct by mapping given cyclicFvsPatchField onto a new patch
		cyclicFvsPatchField
		(
			const cyclicFvsPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, surfaceMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		cyclicFvsPatchField
		(
			const cyclicFvsPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvsPatchField<Type> > clone() const
		{
			return tmp<fvsPatchField<Type> >
			(
				new cyclicFvsPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		cyclicFvsPatchField
		(
			const cyclicFvsPatchField<Type>&,
			const DimensionedField<Type, surfaceMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvsPatchField<Type> > clone
		(
			const DimensionedField<Type, surfaceMesh>& iF
		) const
		{
			return tmp<fvsPatchField<Type> >
			(
				new cyclicFvsPatchField<Type>(*this, iF)
			);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "cyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
