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

Typedef
	Foam::sphericalTensor2D

Description
	SphericalTensor2D of scalars.

SourceFiles
	sphericalTensor2D.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTensor2D_H
#define sphericalTensor2D_H

#include "SphericalTensor2DTemplate.H"
#include "tensor.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef SphericalTensor2D<scalar> sphericalTensor2D;

// Global spherical tensor definitions

static const sphericalTensor2D I2D(1);
static const sphericalTensor2D oneThirdI2D(1.0/3.0);
static const sphericalTensor2D twoThirdsI2D(2.0/3.0);


//- Data associated with sphericalTensor2D type are contiguous
template<>
inline bool contiguous<sphericalTensor2D>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
