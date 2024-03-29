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
	Foam::fvm

Description
	Calculate the BlockLduSystem (matrix & source) for the gradient
	of the field.
	Intended use: block coupled solvers. i.e. implicit grad(p) in momentum
	equation.

Author
	Vuko Vukcevic, FMENA Zagreb.

SourceFiles
	fvmGrad.C

\*---------------------------------------------------------------------------*/

#ifndef fvmGrad_H
#define fvmGrad_H

#include "volFieldsFwd.H"
#include "BlockLduSystem.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


namespace fvm
{
	template<class Type>
	tmp
	<
		BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
	> grad
	(
		const GeometricField<Type, fvPatchField, volMesh>&,
		const word&
	);

	template<class Type>
	tmp
	<
		BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
	> grad
	(
		const GeometricField<Type, fvPatchField, volMesh>&
	);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fvmGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
