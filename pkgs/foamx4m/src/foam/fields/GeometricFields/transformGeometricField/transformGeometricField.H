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

InClass
	Foam::transformGeometricField

Description
	Spatial transformation functions for FieldFields.

SourceFiles
	transformGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef transformGeometricField_H
#define transformGeometricField_H

#include "transform.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
	GeometricField<Type, PatchField, GeoMesh>&,
	const GeometricField<tensor, PatchField, GeoMesh>&,
	const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
	const GeometricField<tensor, PatchField, GeoMesh>&,
	const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
	const GeometricField<tensor, PatchField, GeoMesh>&,
	const tmp<GeometricField<Type, PatchField, GeoMesh> >&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
	const tmp<GeometricField<tensor, PatchField, GeoMesh> >&,
	const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
	const tmp<GeometricField<tensor, PatchField, GeoMesh> >&,
	const tmp<GeometricField<Type, PatchField, GeoMesh> >&
);


template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
	GeometricField<Type, PatchField, GeoMesh>&,
	const dimensionedTensor&,
	const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
	const dimensionedTensor&,
	const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
	const dimensionedTensor&,
	const tmp<GeometricField<Type, PatchField, GeoMesh> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "transformGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
