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

#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<surfaceInterpolationScheme<Type> > scheme
(
	const surfaceScalarField& faceFlux,
	Istream& streamData
)
{
	return surfaceInterpolationScheme<Type>::New
	(
		faceFlux.mesh(),
		faceFlux,
		streamData
	);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<surfaceInterpolationScheme<Type> > scheme
(
	const surfaceScalarField& faceFlux,
	const word& name
)
{
	return surfaceInterpolationScheme<Type>::New
	(
		faceFlux.mesh(),
		faceFlux,
		faceFlux.mesh().schemesDict().interpolationScheme(name)
	);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<surfaceInterpolationScheme<Type> > scheme
(
	const fvMesh& mesh,
	Istream& streamData
)
{
	return surfaceInterpolationScheme<Type>::New
	(
		mesh,
		streamData
	);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<surfaceInterpolationScheme<Type> > scheme
(
	const fvMesh& mesh,
	const word& name
)
{
	return surfaceInterpolationScheme<Type>::New
	(
		mesh,
		mesh.schemesDict().interpolationScheme(name)
	);
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	const surfaceScalarField& faceFlux,
	Istream& schemeData
)
{
	if (surfaceInterpolation::debug)
	{
		Info<< "interpolate"
			<< "(const GeometricField<Type, fvPatchField, volMesh>&, "
			<< "const surfaceScalarField&, Istream&) : "
			<< "interpolating GeometricField<Type, fvPatchField, volMesh> "
			<< endl;
	}

	return scheme<Type>(faceFlux, schemeData)().interpolate(vf);
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	const surfaceScalarField& faceFlux,
	const word& name
)
{
	if (surfaceInterpolation::debug)
	{
		Info<< "interpolate"
			<< "(const GeometricField<Type, fvPatchField, volMesh>&, "
			<< "const surfaceScalarField&, const word&) : "
			<< "interpolating GeometricField<Type, fvPatchField, volMesh> "
			<< "using " << name
			<< endl;
	}

	return scheme<Type>(faceFlux, name)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
	const surfaceScalarField& faceFlux,
	const word& name
)
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf =
		interpolate(tvf(), faceFlux, name);

	tvf.clear();

	return tsf;
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	const tmp<surfaceScalarField>& tFaceFlux,
	const word& name
)
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf =
		interpolate(vf, tFaceFlux(), name);

	tFaceFlux.clear();

	return tsf;
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
	const tmp<surfaceScalarField>& tFaceFlux,
	const word& name
)
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf =
		interpolate(tvf(), tFaceFlux(), name);

	tvf.clear();
	tFaceFlux.clear();

	return tsf;
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	Istream& schemeData
)
{
	if (surfaceInterpolation::debug)
	{
		Info<< "interpolate"
			<< "(const GeometricField<Type, fvPatchField, volMesh>&, "
			<< "Istream&) : "
			<< "interpolating GeometricField<Type, fvPatchField, volMesh> "
			<< endl;
	}

	return scheme<Type>(vf.mesh(), schemeData)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf,
	const word& name
)
{
	if (surfaceInterpolation::debug)
	{
		Info<< "interpolate"
			<< "(const GeometricField<Type, fvPatchField, volMesh>&, "
			<< "const word&) : "
			<< "interpolating GeometricField<Type, fvPatchField, volMesh> "
			<< "using " << name
			<< endl;
	}

	return scheme<Type>(vf.mesh(), name)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
	const word& name
)
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf =
		interpolate(tvf(), name);

	tvf.clear();

	return tsf;
}


// Interpolate field onto faces using central differencing
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	if (surfaceInterpolation::debug)
	{
		Info<< "interpolate"
			<< "(const GeometricField<Type, fvPatchField, volMesh>&) : "
			<< "interpolating GeometricField<Type, fvPatchField, volMesh> "
			<< "using run-time selected scheme"
			<< endl;
	}

	return interpolate(vf, "interpolate(" + vf.name() + ')');
}


// Interpolate field onto faces using central differencing
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
interpolate
(
	const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf =
		interpolate(tvf());
	tvf.clear();
	return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
