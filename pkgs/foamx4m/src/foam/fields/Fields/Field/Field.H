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
	Foam::Field

Description
	Generic templated field type.

SourceFiles
	FieldFunctions.H
	FieldFunctionsM.H
	FieldMapper.H
	FieldM.H
	Field.C
	FieldFunctions.C
	FieldFunctionsM.C

\*---------------------------------------------------------------------------*/

#ifndef Field_H
#define Field_H

#include "tmp.H"
#include "direction.H"
#include "VectorSpace.H"
#include "scalarList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class Field;

//- Pre-declare related SubField type
template<class Type>
class SubField;

template<class Type>
Ostream& operator<<(Ostream&, const Field<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const tmp<Field<Type> >&);

class FieldMapper;
class dictionary;


template<class Type>
class Field
:
	public refCount,
	public List<Type>
{

public:

	// Public typedefs

		//- Field element type (eg. scalar, vector)
		typedef Type PrimitiveType;

		//- Component type
		typedef typename pTraits<Type>::cmptType cmptType;

		//- Declare type of subField
		typedef SubField<Type> subField;


	// Static data members

		static const char* const typeName;

		//- Empty field
		static const Field<Type> zero;


	// Static Member Functions

		//- Return a null field
		inline static const Field<Type>& null()
		{
			return zero;
		}


	// Constructors

		//- Construct null
		//  Used for temporary fields which are initialised after construction
		Field();

		//- Construct given size
		//  Used for temporary fields which are initialised after construction
		explicit Field(const label);

		//- Construct given size and initial value
		Field(const label, const Type&);

		//- Construct as copy of a UList\<Type\>
		explicit Field(const UList<Type>&);

		//- Construct by transferring the List contents
		explicit Field(const Xfer<List<Type> >&);

		//- Construct by 1 to 1 mapping from the given field
		Field
		(
			const UList<Type>& mapF,
			const unallocLabelList& mapAddressing
		);

		//- Construct by 1 to 1 mapping from the given tmp field
		Field
		(
			const tmp<Field<Type> >& tmapF,
			const unallocLabelList& mapAddressing
		);

		//- Construct by interpolative mapping from the given field
		Field
		(
			const UList<Type>& mapF,
			const labelListList& mapAddressing,
			const scalarListList& weights
		);

		//- Construct by interpolative mapping from the given tmp field
		Field
		(
			const tmp<Field<Type> >& tmapF,
			const labelListList& mapAddressing,
			const scalarListList& weights
		);

		//- Construct by mapping from the given field
		Field
		(
			const UList<Type>& mapF,
			const FieldMapper& map
		);

		//- Construct by mapping from the given tmp field
		Field
		(
			const tmp<Field<Type> >& tmapF,
			const FieldMapper& map
		);

		//- Construct as copy
		Field(const Field<Type>&);

		//- Construct as copy or re-use as specified.
		Field(Field<Type>&, bool reUse);

		//- Construct by transferring the Field contents
		Field(const Xfer<Field<Type> >&);

		//- Construct as copy of subField
		Field(const typename Field<Type>::subField&);

		//- Construct as copy of tmp<Field>
#		ifdef ConstructFromTmp
		Field(const tmp<Field<Type> >&);
#		endif

		//- Construct from Istream
		Field(Istream&);

		//- Construct from a dictionary entry
		Field(const word& keyword, const dictionary& dict, const label size);

		//- Clone
		tmp<Field<Type> > clone() const;

		//- Return a pointer to a new Field created on freestore
		static autoPtr<Field<Type> > New(Istream& is)
		{
			return autoPtr<Field<Type> >(new Field<Type>(is));
		}

		//- Return a pointer to a new calculatedFvPatchFieldField created on
		//  freestore without setting patchField values
		template<class Type2>
		static tmp<Field<Type> > NewCalculatedType(const Field<Type2>& f)
		{
			return tmp<Field<Type> >(new Field<Type>(f.size()));
		}


	// Member Functions

		//- 1 to 1 map from the given field
		void map
		(
			const UList<Type>& mapF,
			const unallocLabelList& mapAddressing
		);

		//- 1 to 1 map from the given tmp field
		void map
		(
			const tmp<Field<Type> >& tmapF,
			const unallocLabelList& mapAddressing
		);

		//- Interpolative map from the given field
		void map
		(
			const UList<Type>& mapF,
			const labelListList& mapAddressing,
			const scalarListList& weights
		);

		//- Interpolative map from the given tmp field
		void map
		(
			const tmp<Field<Type> >& tmapF,
			const labelListList& mapAddressing,
			const scalarListList& weights
		);

		//- Map from the given field
		void map
		(
			const UList<Type>& mapF,
			const FieldMapper& map
		);

		//- Map from the given tmp field
		void map
		(
			const tmp<Field<Type> >& tmapF,
			const FieldMapper& map
		);

		//- Map from self
		void autoMap
		(
			const FieldMapper& map
		);

		//- 1 to 1 reverse-map from the given field
		void rmap
		(
			const UList<Type>& mapF,
			const unallocLabelList& mapAddressing
		);

		//- 1 to 1 reverse-map from the given tmp field
		void rmap
		(
			const tmp<Field<Type> >& tmapF,
			const unallocLabelList& mapAddressing
		);

		//- Interpolative reverse map from the given field
		void rmap
		(
			const UList<Type>& mapF,
			const unallocLabelList& mapAddressing,
			const scalarList& weights
		);

		//- Interpolative reverse map from the given tmp field
		void rmap
		(
			const tmp<Field<Type> >& tmapF,
			const unallocLabelList& mapAddressing,
			const scalarList& weights
		);

		//- Negate this field
		void negate();

		//- Return a component field of the field
		tmp<Field<cmptType> > component(const direction) const;

		//- Replace a component field of the field
		void replace(const direction, const UList<cmptType>&);

		//- Replace a component field of the field
		void replace(const direction, const tmp<Field<cmptType> >&);

		//- Replace a component field of the field
		void replace(const direction, const cmptType&);

		//- Return the field transpose (only defined for second rank tensors)
		tmp<Field<Type> > T() const;

		//- Write the field as a dictionary entry
		void writeEntry(const word& keyword, Ostream& os) const;


	// Member operators

		void operator=(const Field<Type>&);
		void operator=(const UList<Type>&);
		void operator=(const SubField<Type>&);
		void operator=(const tmp<Field<Type> >&);
		void operator=(const Type&);

		template<class Form, class Cmpt, int nCmpt>
		void operator=(const VectorSpace<Form,Cmpt,nCmpt>&);

		void operator+=(const UList<Type>&);
		void operator+=(const tmp<Field<Type> >&);

		void operator-=(const UList<Type>&);
		void operator-=(const tmp<Field<Type> >&);

		void operator*=(const UList<scalar>&);
		void operator*=(const tmp<Field<scalar> >&);

		void operator/=(const UList<scalar>&);
		void operator/=(const tmp<Field<scalar> >&);

		void operator+=(const Type&);
		void operator-=(const Type&);

		void operator*=(const scalar&);
		void operator/=(const scalar&);


	// IOstream operators

		friend Ostream& operator<< <Type>
		(Ostream&, const Field<Type>&);

		friend Ostream& operator<< <Type>
		(Ostream&, const tmp<Field<Type> >&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFunctions.H"

#ifdef NoRepository
#	include "Field.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
