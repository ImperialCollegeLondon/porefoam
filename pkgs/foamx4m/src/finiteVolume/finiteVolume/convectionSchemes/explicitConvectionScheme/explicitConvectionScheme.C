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

#include "explicitConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
explicitConvectionScheme<Type>::interpolate
(
	const surfaceScalarField&,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	return tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
explicitConvectionScheme<Type>::flux
(
	const surfaceScalarField& faceFlux,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
explicitConvectionScheme<Type>::fvmDiv
(
	const surfaceScalarField& faceFlux,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	tmp<fvMatrix<Type> > tfvm
	(
		new fvMatrix<Type>
		(
			vf,
			faceFlux.dimensions()*vf.dimensions()
		)
	);
	fvMatrix<Type>& fvm = tfvm();

	// Matrix consistency
	fvm.diag() = 0;

	fvm += this->fvcDiv(faceFlux, vf);

	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfaceFluxCorrection
		= faceFlux*interpolate(faceFlux, vf);

	const fvMesh& mesh = this->mesh();

	if (mesh.schemesDict().fluxRequired(vf.name()))
	{
		fvm.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
	}

	return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
explicitConvectionScheme<Type>::fvcDiv
(
	const surfaceScalarField& faceFlux,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	tmp<GeometricField<Type, fvPatchField, volMesh> > tConvection
	(
		fvc::surfaceIntegrate(flux(faceFlux, vf))
	);

	tConvection().rename
	(
		"convection(" + faceFlux.name() + ',' + vf.name() + ')'
	);

	return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
