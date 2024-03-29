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
	Foam::SlicedGeometricField

Description
	Specialization of GeometricField which holds slices of given complete
	fields in a form that they act as a GeometricField.

	The destructor is wrapped to avoid deallocation of the storage of the
	complete fields when this is destroyed.

	SlicedGeometricField can only be instantiated with a valid form of
	SlicedPatchField to handle the slicing and storage deallocation of the
	boundary field.

SourceFiles
	SlicedGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef SlicedGeometricField_H
#define SlicedGeometricField_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
class SlicedGeometricField
:
	public GeometricField<Type, PatchField, GeoMesh>
{
public:

	typedef typename GeoMesh::Mesh Mesh;
	typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

	class Internal;


private:

	// Private Member Functions

		//- Slice the given field and a create a PtrList of SlicedPatchField
		//  from which the boundary field is built
		tmp<FieldField<PatchField, Type> >  slicedBoundaryField
		(
			const Mesh& mesh,
			const Field<Type>& completeField,
			const bool preserveCouples
		);

		//- Slice the given field and a create a PtrList of SlicedPatchField
		//  from which the boundary field is built
		tmp<FieldField<PatchField, Type> >  slicedBoundaryField
		(
			const Mesh& mesh,
			const FieldField<PatchField, Type>& bField,
			const bool preserveCouples
		);

		////- Disallow default bitwise copy construct
		//SlicedGeometricField(const SlicedGeometricField&);

		//- Disallow default bitwise assignment
		void operator=(const SlicedGeometricField&);

		//- Disallow standard assignment to GeometricField,
		//  == assignment is allowed.
		void operator=(const GeometricField<Type, PatchField, GeoMesh>&);

		//- Disallow standard assignment to tmp<GeometricField>,
		//  == assignment is allowed.
		void operator=(const tmp<GeometricField<Type, PatchField, GeoMesh> >&);


public:

	// Constructors

		//- Construct from components and field to slice
		SlicedGeometricField
		(
			const IOobject&,
			const Mesh&,
			const dimensionSet&,
			const Field<Type>& completeField,
			const bool preserveCouples=true
		);

		//- Construct from components and separate fields to slice for the
		//  internal field and boundary field
		SlicedGeometricField
		(
			const IOobject&,
			const Mesh&,
			const dimensionSet&,
			const Field<Type>& completeIField,
			const Field<Type>& completeBField,
			const bool preserveCouples=true
		);

		//- Construct from GeometricField. Reuses full internal and
		//  patch fields except on couples (preserveCouples=true).
		SlicedGeometricField
		(
			const IOobject&,
			const GeometricField<Type, PatchField, GeoMesh>&,
			const bool preserveCouples=true
		);

		//- Construct as copy
		SlicedGeometricField
		(
			const SlicedGeometricField
			<
				Type,
				PatchField,
				SlicedPatchField,
				GeoMesh
			>&
		);


	//- Destructor
	~SlicedGeometricField();


	// Member functions

		//- Reset storage of a field to slice
		void reset(const Field<Type>& completeField);

		//- Reset storage of a field to slice with separate fields for the
		//  internal field and boundary field
		void reset
		(
			const Field<Type>& completeIField,
			const Field<Type>& completeBField
		);

		//- Correct boundary field
		void correctBoundaryConditions();
};



//- The internalField of a SlicedGeometricField
template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
class SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
Internal
:
	public GeometricField<Type, PatchField, GeoMesh>::Internal
{

public:

	// Constructors

		//- Construct from components and field to slice
		Internal
		(
			const IOobject&,
			const Mesh&,
			const dimensionSet&,
			const Field<Type>& iField
		);


	//- Destructor
	~Internal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "SlicedGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
