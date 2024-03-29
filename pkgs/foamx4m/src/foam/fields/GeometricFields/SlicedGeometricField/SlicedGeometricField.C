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

\*---------------------------------------------------------------------------*/

#include "SlicedGeometricField.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::tmp<Foam::FieldField<PatchField, Type> >
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
slicedBoundaryField
(
	const Mesh& mesh,
	const Field<Type>& completeField,
	const bool preserveCouples
)
{
	tmp<FieldField<PatchField, Type> > tbf
	(
		new FieldField<PatchField, Type>(mesh.boundary().size())
	);

	FieldField<PatchField, Type>& bf = tbf();

	forAll (mesh.boundary(), patchi)
	{
		if (preserveCouples && mesh.boundary()[patchi].coupled())
		{
			// For coupled patches construct the correct patch field type
			// This is a normal patch field where we can assign the values
			// Bug fix: New will already create the correct type on the boundary
			// HJ, 4/Jan/2009
			bf.set
			(
				patchi,
				PatchField<Type>::New
				(
					PatchField<Type>::calculatedType(),
					mesh.boundary()[patchi],
					*this
				)
			);

			// Initialize the values on the coupled patch to those of the slice
			// of the given field.
			// Note: these will usually be over-ridden by the boundary field
			// evaluation e.g. in the case of processor and cyclic patches.
			bf[patchi] = SlicedPatchField<Type>
			(
				mesh.boundary()[patchi],
				DimensionedField<Type, GeoMesh>::null(),
				completeField
			);
		}
		else
		{
			// Create a sliced patch field, referencing into external data
			bf.set
			(
				patchi,
				new SlicedPatchField<Type>
				(
					mesh.boundary()[patchi],
					DimensionedField<Type, GeoMesh>::null(),
					completeField
				)
			);
		}
	}

	return tbf;
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::tmp<Foam::FieldField<PatchField, Type> >
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
slicedBoundaryField
(
	const Mesh& mesh,
	const FieldField<PatchField, Type>& bField,
	const bool preserveCouples
)
{
	tmp<FieldField<PatchField, Type> > tbf
	(
		new FieldField<PatchField, Type>(mesh.boundary().size())
	);

	FieldField<PatchField, Type>& bf = tbf();

	forAll (mesh.boundary(), patchi)
	{
		if (preserveCouples && mesh.boundary()[patchi].coupled())
		{
			// For coupled patches construct the correct patch field type
			bf.set
			(
				patchi,
				PatchField<Type>::New
				(
					mesh.boundary()[patchi].type(),
					mesh.boundary()[patchi],
					*this
				)
			);

			// Assign field
			bf[patchi] == bField[patchi];
		}
		else
		{
			// Create unallocated copy of patch field
			bf.set
			(
				patchi,
				new SlicedPatchField<Type>
				(
					mesh.boundary()[patchi],
					DimensionedField<Type, GeoMesh>::null()
				)
			);
			bf[patchi].UList<Type>::operator=(bField[patchi]);
		}
	}

	return tbf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
Internal::Internal
(
	const IOobject& io,
	const Mesh& mesh,
	const dimensionSet& ds,
	const Field<Type>& iField
)
:
	DimensionedField<Type, GeoMesh>
	(
		io,
		mesh,
		ds,
		Field<Type>()
	)
{
	// Set the internalField to the slice of the complete field
	UList<Type>::operator=
	(
		typename Field<Type>::subField(iField, GeoMesh::size(mesh))
	);
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
	const IOobject& io,
	const Mesh& mesh,
	const dimensionSet& ds,
	const Field<Type>& completeField,
	const bool preserveCouples
)
:
	GeometricField<Type, PatchField, GeoMesh>
	(
		io,
		mesh,
		ds,
		Field<Type>(),
		slicedBoundaryField(mesh, completeField, preserveCouples)
	)
{
	// Set the internalField to the slice of the complete field
	UList<Type>::operator=
	(
		typename Field<Type>::subField(completeField, GeoMesh::size(mesh))
	);

	correctBoundaryConditions();
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
	const IOobject& io,
	const Mesh& mesh,
	const dimensionSet& ds,
	const Field<Type>& completeIField,
	const Field<Type>& completeBField,
	const bool preserveCouples
)
:
	GeometricField<Type, PatchField, GeoMesh>
	(
		io,
		mesh,
		ds,
		Field<Type>(),
		slicedBoundaryField(mesh, completeBField, preserveCouples)
	)
{
	// Set the internalField to the slice of the complete field
	UList<Type>::operator=
	(
		typename Field<Type>::subField(completeIField, GeoMesh::size(mesh))
	);

	correctBoundaryConditions();
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
	const IOobject& io,
	const GeometricField<Type, PatchField, GeoMesh>& gf,
	const bool preserveCouples
)
:
	GeometricField<Type, PatchField, GeoMesh>
	(
		io,
		gf.mesh(),
		gf.dimensions(),
		Field<Type>(),
		slicedBoundaryField(gf.mesh(), gf.boundaryField(), preserveCouples)
	)
{
	// Set the internalField to the supplied internal field
	UList<Type>::operator=(gf.internalField());

	correctBoundaryConditions();
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
SlicedGeometricField
(
	const SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>& gf
)
:
	GeometricField<Type, PatchField, GeoMesh>
	(
		gf,
		gf.mesh(),
		gf.dimensions(),
		Field<Type>(),
		slicedBoundaryField(gf.mesh(), gf.boundaryField(), true)
	)
{
	// Set the internalField to the supplied internal field
	UList<Type>::operator=(gf.internalField());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
~SlicedGeometricField()
{
	// Set the internalField storage pointer to nullptr before its destruction
	// to protect the field it a slice of.
	UList<Type>::operator=(UList<Type>(nullptr, 0));
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
Internal::~Internal()
{
	// Set the internalField storage pointer to nullptr before its destruction
	// to protect the field it a slice of.
	UList<Type>::operator=(UList<Type>(nullptr, 0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
void
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::reset
(
	const Field<Type>& completeField
)
{
	// Set the internalField to the slice of the complete field
	UList<Type>::operator=
	(
		typename Field<Type>::subField(completeField, this->size())
	);

	FieldField<PatchField, Type>& bf = this->boundaryField();

	const fvBoundaryMesh& bMesh = this->mesh().boundary();

	forAll (bMesh, patchi)
	{
		// Note: assuming preserveCouples = true
		// HJ, 1/Dec/2017
		if (bMesh[patchi].coupled())
		{
			// Initialize the values on the coupled patch to those of the slice
			// of the given field.
			// Note: these will usually be over-ridden by the boundary field
			// evaluation e.g. in the case of processor and cyclic patches.
			bf[patchi] = SlicedPatchField<Type>
			(
				bMesh[patchi],
				DimensionedField<Type, GeoMesh>::null(),
				completeField
			);
		}
		else
		{
			bf[patchi].UList<Type>::operator=
			(
				bMesh[patchi].patchSlice(completeField)
			);
		}
	}
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
void
Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::reset
(
	const Field<Type>& completeIField,
	const Field<Type>& completeBField
)
{
	// Set the internalField to the slice of the complete field
	UList<Type>::operator=
	(
		typename Field<Type>::subField(completeIField, this->size())
	);

	FieldField<PatchField, Type>& bf = this->boundaryField();

	const fvBoundaryMesh& bMesh = this->mesh().boundary();

	forAll (bMesh, patchi)
	{
		// Note: assuming preserveCouples = true
		// HJ, 1/Dec/2017
		if (bMesh[patchi].coupled())
		{
			// Assign field
			bf[patchi] = SlicedPatchField<Type>
			(
				bMesh[patchi],
				DimensionedField<Type, GeoMesh>::null(),
				completeBField
			);
		}
		else
		{
			bf[patchi].UList<Type>::operator=
			(
				bMesh[patchi].patchSlice(completeBField)
			);
		}
	}
}


template
<
	class Type,
	template<class> class PatchField,
	template<class> class SlicedPatchField,
	class GeoMesh
>
void Foam::SlicedGeometricField<Type, PatchField, SlicedPatchField, GeoMesh>::
correctBoundaryConditions()
{
	GeometricField<Type, PatchField, GeoMesh>::correctBoundaryConditions();
}


// ************************************************************************* //
