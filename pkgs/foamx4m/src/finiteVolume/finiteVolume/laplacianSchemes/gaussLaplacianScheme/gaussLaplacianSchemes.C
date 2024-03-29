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
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
	makeFvLaplacianScheme(gaussLaplacianScheme)
}
}

#define declareFvmLaplacianScalarGamma(Type)                                 \
					                                                         \
template<>                                                                   \
Foam::tmp<Foam::fvMatrix<Foam::Type> >                                       \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvmLaplacian       \
(                                                                            \
	const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,         \
	const GeometricField<Type, fvPatchField, volMesh>& vf                    \
)                                                                            \
{                                                                            \
	const fvMesh& mesh = this->mesh();                                       \
					                                                         \
	GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf            \
	(                                                                        \
		gamma*mesh.magSf()                                                   \
	);                                                                       \
					                                                         \
	tmp<fvMatrix<Type> > tfvm = fvmLaplacianUncorrected(gammaMagSf, vf);     \
	fvMatrix<Type>& fvm = tfvm();                                            \
					                                                         \
	if (this->tsnGradScheme_().corrected())                                  \
	{                                                                        \
		if (mesh.schemesDict().fluxRequired(vf.name()))                      \
		{                                                                    \
			fvm.faceFluxCorrectionPtr() = new                                \
			GeometricField<Type, fvsPatchField, surfaceMesh>                 \
			(                                                                \
				gammaMagSf*this->tsnGradScheme_().correction(vf)             \
			);                                                               \
					                                                         \
			fvm.source() -=                                                  \
				mesh.V()*                                                    \
				fvc::div                                                     \
				(                                                            \
					*fvm.faceFluxCorrectionPtr()                             \
				)().internalField();                                         \
		}                                                                    \
		else                                                                 \
		{                                                                    \
			fvm.source() -=                                                  \
				mesh.V()*                                                    \
				fvc::div                                                     \
				(                                                            \
					gammaMagSf*this->tsnGradScheme_().correction(vf)         \
				)().internalField();                                         \
		}                                                                    \
	}                                                                        \
					                                                         \
	return tfvm;                                                             \
}                                                                            \
					                                                         \
					                                                         \
template<>                                                                   \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh> >\
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvcLaplacian       \
(                                                                            \
	const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,         \
	const GeometricField<Type, fvPatchField, volMesh>& vf                    \
)                                                                            \
{                                                                            \
	const fvMesh& mesh = this->mesh();                                       \
					                                                         \
	tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian             \
	(                                                                        \
		fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())       \
	);                                                                       \
					                                                         \
	tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');\
					                                                         \
	return tLaplacian;                                                       \
}


declareFvmLaplacianScalarGamma(scalar);
declareFvmLaplacianScalarGamma(vector);
declareFvmLaplacianScalarGamma(sphericalTensor);
declareFvmLaplacianScalarGamma(symmTensor);
declareFvmLaplacianScalarGamma(tensor);


// ************************************************************************* //
