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

Class
	Foam::PhiLimiter

Description
	Class with limiter function which returns the limiter for the
	Phi differencing scheme.

	Used in conjunction with the template class PhiScheme.

SourceFiles
	Phi.C

\*---------------------------------------------------------------------------*/

#ifndef Phi_H
#define Phi_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class PhiLimiter
{
	scalar k_;

public:

	PhiLimiter(Istream& is)
	:
		k_(readScalar(is))
	{
		if (k_ < 0 || k_ > 1)
		{
			FatalIOErrorIn("PhiLimiter(Istream& is)", is)
				<< "coefficient = " << k_
				<< " should be >= 0 and <= 1"
				<< exit(FatalIOError);
		}
	}

	scalar limiter
	(
		const scalar cdWeight,
		const scalar faceFlux,
		const vector& PhiP,
		const vector& PhiN,
		const vector& Sf,
		const scalar&
	) const
	{
		scalar phiP = Sf&PhiP;
		scalar phiN = Sf&PhiN;

		scalar phiU;

		if (faceFlux > 0)
		{
			phiU = phiP;
		}
		else
		{
			phiU = phiN;
		}

		scalar phiCD = cdWeight*phiP + (1 - cdWeight)*phiN;

		// Calculate the effective limiter for the Phi interpolation
		//scalar PLimiter =
		//    (1.0 - k_) + k_*(faceFlux - phiU)/stabilise(phiCD - phiU, SMALL);

		scalar PLimiter =
			((faceFlux - phiU)/stabilise(phiCD - phiU, SMALL) + k_);

		// Limit the limiter between upwind and central
		return max(min(PLimiter, 1), 0);
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
