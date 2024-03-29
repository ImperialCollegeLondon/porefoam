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

Description
	Abstract base class for surface interpolation schemes.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from face-flux field and coefficient
template<class Type>
multivariateSurfaceInterpolationScheme<Type>::
multivariateSurfaceInterpolationScheme
(
	const fvMesh& mesh,
	const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
	const surfaceScalarField&,
	Istream&
)
:
	mesh_(mesh),
	fields_(vtfs)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<multivariateSurfaceInterpolationScheme<Type> >
multivariateSurfaceInterpolationScheme<Type>::New
(
	const fvMesh& mesh,
	const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
	const surfaceScalarField& faceFlux,
	Istream& schemeData
)
{
	if (fv::debug)
	{
		Info<< "multivariateSurfaceInterpolationScheme<Type>::New"
			   "(const fvMesh& mesh, const fieldTable&, "
			   "const surfaceScalarField&, Istream&) : "
			   "constructing surfaceInterpolationScheme<Type>"
			<< endl;
	}

	const word schemeName(schemeData);

	typename IstreamConstructorTable::iterator constructorIter =
		IstreamConstructorTablePtr_->find(schemeName);

	if (constructorIter == IstreamConstructorTablePtr_->end())
	{
		FatalIOErrorIn
		(
			"multivariateSurfaceInterpolationScheme<Type>::New"
			"(const fvMesh& mesh, const fieldTable&, "
			"const surfaceScalarField&, Istream&)",
			schemeData
		)   << "unknown discretisation scheme " << schemeName << endl << endl
			<< "Valid schemes are :" << endl
			<< IstreamConstructorTablePtr_->sortedToc()
			<< exit(FatalIOError);
	}

	return constructorIter()(mesh, vtfs, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
multivariateSurfaceInterpolationScheme<Type>::
~multivariateSurfaceInterpolationScheme()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
