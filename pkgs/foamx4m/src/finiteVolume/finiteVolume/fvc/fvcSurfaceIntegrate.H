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

InNamespace
	Foam::fvc

Description
	Surface integrate surfaceField creating a volField.
	Surface sum a surfaceField creating a volField.

SourceFiles
	fvcSurfaceIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fvcSurfaceIntegrate_H
#define fvcSurfaceIntegrate_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


namespace fvc
{
	template<class Type>
	void surfaceIntegrate
	(
		Field<Type>&,
		const GeometricField<Type, fvsPatchField, surfaceMesh>&
	);

	template<class Type>
	tmp<GeometricField<Type, fvPatchField, volMesh> >
	surfaceIntegrate
	(
		const GeometricField<Type, fvsPatchField, surfaceMesh>&
	);

	template<class Type>
	tmp<GeometricField<Type, fvPatchField, volMesh> >
	surfaceIntegrate
	(
		const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
	);

	template<class Type>
	tmp<GeometricField<Type, fvPatchField, volMesh> > surfaceSum
	(
		const GeometricField<Type, fvsPatchField, surfaceMesh>&
	);

	template<class Type>
	tmp<GeometricField<Type, fvPatchField, volMesh> > surfaceSum
	(
		const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
	);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fvcSurfaceIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
