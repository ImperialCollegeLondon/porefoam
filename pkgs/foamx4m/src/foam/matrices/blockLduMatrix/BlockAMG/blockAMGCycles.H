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

#include "BlockAMGCycle.H"

#ifdef USE_VecN
// KRJ included to make the different types
#include "VectorTensorNFieldsFwd.H"
#include "ExpandTensorNField.H"
#include "VectorNFieldTypes.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// KRJ: 2012-12-15: Is this enough to have the templating?

typedef BlockAMGCycle<scalar> blockAMGScalarCycle;
typedef BlockAMGCycle<vector> blockAMGVectorCycle;
typedef BlockAMGCycle<tensor> blockAMGTensorCycle;

#ifdef USE_VecN
#define makeTypedef(type, Type, args...)                \
	typedef BlockAMGCycle<type > blockAMG##Type##Cycle;

forAllVectorNTypes(makeTypedef)
#undef makeTypedef
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
