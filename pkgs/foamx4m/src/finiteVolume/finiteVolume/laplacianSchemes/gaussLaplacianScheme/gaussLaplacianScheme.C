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

#include "gaussLaplacianScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type> >
gaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
(
	const surfaceScalarField& gammaMagSf,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	tmp<surfaceScalarField> tdeltaCoeffs =
		this->tsnGradScheme_().deltaCoeffs(vf);
	const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

	tmp<fvMatrix<Type> > tfvm
	(
		new fvMatrix<Type>
		(
			vf,
			deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
		)
	);
	fvMatrix<Type>& fvm = tfvm();

	fvm.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField();
	fvm.negSumDiag();

	forAll(fvm.psi().boundaryField(), patchI)
	{
		const fvPatchField<Type>& psf = fvm.psi().boundaryField()[patchI];
		const fvsPatchScalarField& patchGamma =
			gammaMagSf.boundaryField()[patchI];

		fvm.internalCoeffs()[patchI] = patchGamma*psf.gradientInternalCoeffs();
		fvm.boundaryCoeffs()[patchI] = -patchGamma*psf.gradientBoundaryCoeffs();
	}

	// Manipulate internal and boundary coeffs for diffusion. Needed for very
	// special treatment and is currently used only for ensuring implicit
	// conservation across GGI interface that has partially covered faces. Does
	// nothing for other fvPatchFields. VV, 8/Mar/2018.
	forAll(fvm.psi().boundaryField(), patchI)
	{
		fvm.psi().boundaryField()[patchI].manipulateGradientCoeffs(fvm);
	}

	return tfvm;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
gaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
	const surfaceVectorField& SfGammaCorr,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	const fvMesh& mesh = this->mesh();

	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tgammaSnGradCorr
	(
		new GeometricField<Type, fvsPatchField, surfaceMesh>
		(
			IOobject
			(
				"gammaSnGradCorr("+vf.name()+')',
				vf.instance(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			SfGammaCorr.dimensions()*
			vf.dimensions()*mesh.deltaCoeffs().dimensions()
		)
	);

	for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
	{
		tgammaSnGradCorr().replace
		(
			cmpt,
			SfGammaCorr & fvc::interpolate(fvc::grad(vf.component(cmpt)))
		);
	}

	return tgammaSnGradCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussLaplacianScheme<Type, GType>::fvcLaplacian
(
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	const fvMesh& mesh = this->mesh();

	tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
	(
		fvc::div(this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
	);

	tLaplacian().rename("laplacian(" + vf.name() + ')');

	return tLaplacian;
}


template<class Type, class GType>
tmp<fvMatrix<Type> >
gaussLaplacianScheme<Type, GType>::fvmLaplacian
(
	const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	const fvMesh& mesh = this->mesh();

	const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());

	const surfaceVectorField SfGamma(mesh.Sf() & gamma);
	const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
	(
		SfGamma & Sn
	);
	const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

	tmp<fvMatrix<Type> > tfvm = fvmLaplacianUncorrected(SfGammaSn, vf);
	fvMatrix<Type>& fvm = tfvm();

	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfaceFluxCorrection
		= gammaSnGradCorr(SfGammaCorr, vf);

	if (this->tsnGradScheme_().corrected())
	{
		tfaceFluxCorrection() +=
			SfGammaSn*this->tsnGradScheme_().correction(vf);
	}

	fvm.source() -= mesh.V()*fvc::div(tfaceFluxCorrection())().internalField();

	if (mesh.schemesDict().fluxRequired(vf.name()))
	{
		fvm.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
	}

	return tfvm;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussLaplacianScheme<Type, GType>::fvcLaplacian
(
	const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	const fvMesh& mesh = this->mesh();

	const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
	const surfaceVectorField SfGamma(mesh.Sf() & gamma);
	const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
	(
		SfGamma & Sn
	);
	const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

	tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
	(
		fvc::div
		(
			SfGammaSn*this->tsnGradScheme_().snGrad(vf)
		  + gammaSnGradCorr(SfGammaCorr, vf)
		)
	);

	tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');

	return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
