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

#include "gaussConvectionScheme.H"
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
gaussConvectionScheme<Type>::interpolate
(
	const surfaceScalarField&,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	return tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
gaussConvectionScheme<Type>::flux
(
	const surfaceScalarField& faceFlux,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
gaussConvectionScheme<Type>::fvmDiv
(
	const surfaceScalarField& faceFlux,
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
	const surfaceScalarField& weights = tweights();

	tmp<fvMatrix<Type> > tfvm
	(
		new fvMatrix<Type>
		(
			vf,
			faceFlux.dimensions()*vf.dimensions()
		)
	);
	fvMatrix<Type>& fvm = tfvm();

	fvm.lower() = -weights.internalField()*faceFlux.internalField();
	fvm.upper() = fvm.lower() + faceFlux.internalField();
	fvm.negSumDiag();

	forAll(fvm.psi().boundaryField(), patchI)
	{
		const fvPatchField<Type>& psf = fvm.psi().boundaryField()[patchI];
		const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchI];
		const fvsPatchScalarField& pw = weights.boundaryField()[patchI];

		fvm.internalCoeffs()[patchI] = patchFlux*psf.valueInternalCoeffs(pw);
		fvm.boundaryCoeffs()[patchI] = -patchFlux*psf.valueBoundaryCoeffs(pw);
	}

	// Manipulate internal and boundary coeffs for convection. Needed for very
	// special treatment and is currently used only for ensuring implicit
	// conservation across GGI interface that has partially covered faces. Does
	// nothing for other fvPatchFields. VV, 8/Mar/2018.
	forAll(fvm.psi().boundaryField(), patchI)
	{
		fvm.psi().boundaryField()[patchI].manipulateValueCoeffs(fvm);
	}

	if (tinterpScheme_().corrected())
	{
		fvm += fvc::surfaceIntegrate(faceFlux*tinterpScheme_().correction(vf));
	}

	return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussConvectionScheme<Type>::fvcDiv
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
