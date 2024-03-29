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

#include "volPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "globalPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void volPointInterpolation::interpolateInternalField
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
	if (debug)
	{
		Info<< "volPointInterpolation::interpolateInternalField("
			<< "const GeometricField<Type, fvPatchField, volMesh>&, "
			<< "GeometricField<Type, pointPatchField, pointMesh>&) : "
			<< "interpolating field from cells to points"
			<< endl;
	}

	const labelListList& pointCells = vf.mesh().pointCells();

	// Get point weights
	const scalarListList& pWeights = this->pointWeights();

	// Multiply volField by weighting factor matrix to create pointField
	forAll(pointCells, pointi)
	{
		const scalarList& pw = pWeights[pointi];
		const labelList& ppc = pointCells[pointi];

		pf[pointi] = pTraits<Type>::zero;

		forAll(ppc, pointCelli)
		{
			pf[pointi] += pw[pointCelli]*vf[ppc[pointCelli]];
		}
	}
}


template<class Type>
void volPointInterpolation::interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
	if (debug)
	{
		Info<< "volPointInterpolation::interpolate("
			<< "const GeometricField<Type, fvPatchField, volMesh>&, "
			<< "GeometricField<Type, pointPatchField, pointMesh>&) : "
			<< "interpolating field from cells to points"
			<< endl;
	}

	interpolateInternalField(vf, pf);

	// Interpolate to the patches preserving fixed value BCs
	boundaryInterpolator_.interpolate(vf, pf, false);
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	const wordList& patchFieldTypes
) const
{
	wordList types(patchFieldTypes);

	const pointMesh& pMesh = pointMesh::New(vf.mesh());

	// If the last patch of the pointBoundaryMesh is the global patch
	// it must be added to the list of patchField types
	if
	(
		isType<globalPointPatch>
		(
			pMesh.boundary()[pMesh.boundary().size() - 1]
		)
	)
	{
		types.setSize(types.size() + 1);
		types[types.size()-1] = pMesh.boundary()[types.size()-1].type();
	}

	// Construct tmp<pointField>
	tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
	(
		new GeometricField<Type, pointPatchField, pointMesh>
		(
			IOobject
			(
				"volPointInterpolate(" + vf.name() + ')',
				vf.instance(),
				pMesh.thisDb()
			),
			pMesh,
			vf.dimensions(),
			types
		)
	);

	interpolateInternalField(vf, tpf());

	// Interpolate to the patches overriding fixed value BCs
	boundaryInterpolator_.interpolate(vf, tpf(), true);

	return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
	const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
	const wordList& patchFieldTypes
) const
{
	// Construct tmp<pointField>
	tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
		interpolate(tvf(), patchFieldTypes);
	tvf.clear();
	return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	const pointMesh& pm = pointMesh::New(vf.mesh());

	tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
	(
		new GeometricField<Type, pointPatchField, pointMesh>
		(
			IOobject
			(
				"volPointInterpolate(" + vf.name() + ')',
				vf.instance(),
				pm.thisDb()
			),
			pm,
			vf.dimensions()
		)
	);

	interpolateInternalField(vf, tpf());
	boundaryInterpolator_.interpolate(vf, tpf(), false);

	return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
	const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
) const
{
	// Construct tmp<pointField>
	tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
		interpolate(tvf());
	tvf.clear();
	return tpf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
