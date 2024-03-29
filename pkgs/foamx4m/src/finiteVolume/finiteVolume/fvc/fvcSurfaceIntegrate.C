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

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
	Field<Type>& ivf,
	const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
	const fvMesh& mesh = ssf.mesh();

	const unallocLabelList& owner = mesh.owner();
	const unallocLabelList& neighbour = mesh.neighbour();

	const Field<Type>& issf = ssf;

	forAll(owner, facei)
	{
		ivf[owner[facei]] += issf[facei];
		ivf[neighbour[facei]] -= issf[facei];
	}

	forAll(mesh.boundary(), patchi)
	{
		const unallocLabelList& pFaceCells =
			mesh.boundary()[patchi].faceCells();

		const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

		forAll(mesh.boundary()[patchi], facei)
		{
			ivf[pFaceCells[facei]] += pssf[facei];
		}
	}

	ivf /= mesh.V();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
	const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
	const fvMesh& mesh = ssf.mesh();

	tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
	(
		new GeometricField<Type, fvPatchField, volMesh>
		(
			IOobject
			(
				"surfaceIntegrate("+ssf.name()+')',
				ssf.instance(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensioned<Type>
			(
				"0",
				ssf.dimensions()/dimVol,
				pTraits<Type>::zero
			),
			zeroGradientFvPatchField<Type>::typeName
		)
	);
	GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

	surfaceIntegrate(vf.internalField(), ssf);
	vf.correctBoundaryConditions();

	return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
	const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
	tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
	(
		fvc::surfaceIntegrate(tssf())
	);
	tssf.clear();
	return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceSum
(
	const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
	const fvMesh& mesh = ssf.mesh();

	tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
	(
		new GeometricField<Type, fvPatchField, volMesh>
		(
			IOobject
			(
				"surfaceSum("+ssf.name()+')',
				ssf.instance(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero),
			zeroGradientFvPatchField<Type>::typeName
		)
	);
	GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

	const unallocLabelList& owner = mesh.owner();
	const unallocLabelList& neighbour = mesh.neighbour();

	forAll(owner, facei)
	{
		vf[owner[facei]] += ssf[facei];
		vf[neighbour[facei]] += ssf[facei];
	}

	forAll(mesh.boundary(), patchi)
	{
		const unallocLabelList& pFaceCells =
			mesh.boundary()[patchi].faceCells();

		const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

		forAll(mesh.boundary()[patchi], facei)
		{
			vf[pFaceCells[facei]] += pssf[facei];
		}
	}

	vf.correctBoundaryConditions();

	return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > surfaceSum
(
	const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
	tmp<GeometricField<Type, fvPatchField, volMesh> > tvf = surfaceSum(tssf());
	tssf.clear();
	return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
