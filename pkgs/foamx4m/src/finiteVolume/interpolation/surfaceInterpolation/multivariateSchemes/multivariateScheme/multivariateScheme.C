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

#include "volFields.H"
#include "surfaceFields.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Scheme>
multivariateScheme<Type, Scheme>::multivariateScheme
(
	const fvMesh& mesh,
	const typename multivariateSurfaceInterpolationScheme<Type>::
		fieldTable& fields,
	const surfaceScalarField& faceFlux,
	Istream& schemeData
)
:
	multivariateSurfaceInterpolationScheme<Type>
	(
		mesh,
		fields,
		faceFlux,
		schemeData
	),
	Scheme::LimiterType(schemeData),
	faceFlux_(faceFlux),
	weights_
	(
		IOobject
		(
			"multivariateWeights",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimless
	)
{
	typename multivariateSurfaceInterpolationScheme<Type>::
		fieldTable::const_iterator iter = this->fields().begin();

	surfaceScalarField limiter =
		Scheme(mesh, faceFlux_, *this).limiter(*iter());

	for (++iter; iter != this->fields().end(); ++iter)
	{
		limiter = min
		(
			limiter,
			Scheme(mesh, faceFlux_, *this).limiter(*iter())
		);
	}

	weights_ =
		limiter*mesh.surfaceInterpolation::weights()
	  + (scalar(1) - limiter)*upwind<Type>(mesh, faceFlux_).weights();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
