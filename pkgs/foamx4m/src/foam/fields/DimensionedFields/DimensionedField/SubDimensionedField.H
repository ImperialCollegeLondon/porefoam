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
	Foam::SubDimensionedField

Description
	SubDimensionedField is a DimensionedField obtained as a section of another
	DimensionedField.

	Thus it is itself unallocated so that no storage is allocated or
	deallocated during its use.  To achieve this behaviour,
	SubDimensionedField is derived from SubField rather than Field.

SourceFiles
	SubDimensionedFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef SubDimensionedField_H
#define SubDimensionedField_H

#include "Field.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type, class GeoMesh>
class SubDimensionedField
:
	public regIOobject,
	public SubField<Type>
{

public:

	// Public typedefs

		typedef typename GeoMesh::Mesh Mesh;
		typedef typename Field<Type>::cmptType cmptType;


	// Static Member Functions

		//- Return a null SubDimensionedField
		//  Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
		static inline const SubDimensionedField<Type, GeoMesh>& null();


	// Constructors

		//- Construct from a SubField
		inline SubDimensionedField
		(
			const SubField<Type>& slist
		);

		//- Construct from a UList and size
		inline SubDimensionedField
		(
			const UList<Type>& list,
			const label subSize
		);

		//- Construct from a UList start and end indices
		inline SubDimensionedField
		(
			const UList<Type>& list,
			const label subSize,
			const label startIndex
		);

		//- Construct as copy
		inline SubDimensionedField
		(
			const SubDimensionedField<cmptType, GeoMesh>& sfield
		);


	// Member functions

		//- Return a component field of the field
		inline tmp<DimensionedField<cmptType, GeoMesh> > component
		(
			const direction
		) const;

		//- Return the field transpose (only defined for second rank tensors)
		tmp<DimensionedField<Type, GeoMesh> > T() const;


	// Member operators

		//- Assignment
		inline void operator=(const SubDimensionedField<Type, GeoMesh>&);

		//- Allow cast to a const DimensionedField<Type, GeoMesh>&
		inline operator const DimensionedField<Type, GeoMesh>&() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubDimensionedFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
