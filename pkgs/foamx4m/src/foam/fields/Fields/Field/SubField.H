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
	Foam::SubField

Description
	SubField is a Field obtained as a section of another Field.

	Thus it is itself unallocated so that no storage is allocated or
	deallocated during its use.  To achieve this behaviour, SubField is
	derived from a SubList rather than a List.

SourceFiles
	SubFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef SubField_H
#define SubField_H

#include "SubList.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Pre-declare SubField and related Field type
template<class Type> class Field;
template<class Type> class SubField;


template<class Type>
class SubField
:
	public refCount,
	public SubList<Type>
{

public:

	//- Component type
	typedef typename pTraits<Type>::cmptType cmptType;


	// Static Member Functions

		//- Return a null SubField
		//  Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
		static inline const SubField<Type>& null();


	// Constructors

		//- Construct from a SubList
		inline SubField(const SubList<Type>&);

		//- Construct from a UList\<Type\>, using the entire size
		explicit inline SubField(const UList<Type>&);

		//- Construct from a UList\<Type\> with a given size
		inline SubField
		(
			const UList<Type>& list,
			const label subSize
		);

		//- Construct from a UList\<Type\> with a given size and start index
		inline SubField
		(
			const UList<Type>& list,
			const label subSize,
			const label startIndex
		);

		//- Construct as copy
		inline SubField(const SubField<Type>&);


	// Member functions

		//- Return a component field of the field
		inline tmp<Field<cmptType> > component(const direction) const;

		//- Return the field transpose (only defined for second rank tensors)
		tmp<Field<Type> > T() const;


	// Member operators

		//- Assignment via UList operator. Takes linear time.
		inline void operator=(const SubField<Type>&);

		//- Allow cast to a const Field\<Type\>&
		inline operator const Field<Type>&() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
