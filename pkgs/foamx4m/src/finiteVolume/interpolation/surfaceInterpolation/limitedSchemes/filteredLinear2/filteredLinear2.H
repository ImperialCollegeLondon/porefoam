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
	Foam::filteredLinear2Limiter

Description
	Class to generate weighting factors for the filtered-linear-2
	differencing scheme.

	The aim is to remove high-frequency modes with "staggering"
	characteristics by comparing the face gradient with both neighbouring
	cell gradients and introduce small amounts of upwind in order to damp
	these modes.

	Used in conjunction with the template class LimitedScheme.

SourceFiles
	filteredLinear2.C

\*---------------------------------------------------------------------------*/

#ifndef filteredLinear2_H
#define filteredLinear2_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class LimiterFunc>
class filteredLinear2Limiter
:
	public LimiterFunc
{
	// Private data

		// Scaling corefficient for the gradient ratio,
		// 0 = linear
		// 1 = fully limited
		scalar k_;

		// Maximum allowed overshoot/undershoot relative to the difference
		// across the face.
		// On input:
		//     0 = no overshoot/undershoot
		//     1 = overshoot/undershoot equal to the difference across the face
		// Note: After input 1 is added to l_
		scalar l_;


public:

	filteredLinear2Limiter(Istream& is)
	:
		k_(readScalar(is)),
		l_(readScalar(is))
	{
		if (k_ < 0 || k_ > 1)
		{
			FatalIOErrorIn("filteredLinearV2Limiter(Istream& is)", is)
				<< "coefficient = " << k_
				<< " should be >= 0 and <= 1"
				<< exit(FatalIOError);
		}

		if (l_ < 0 || l_ > 1)
		{
			FatalIOErrorIn("filteredLinearV2Limiter(Istream& is)", is)
				<< "coefficient = " << l_
				<< " should be >= 0 and <= 1"
				<< exit(FatalIOError);
		}

		l_ += 1.0;
	}

	scalar limiter
	(
		const scalar cdWeight,
		const scalar faceFlux,
		const typename LimiterFunc::phiType& phiP,
		const typename LimiterFunc::phiType& phiN,
		const typename LimiterFunc::gradPhiType& gradcP,
		const typename LimiterFunc::gradPhiType& gradcN,
		const vector& d
	) const
	{
		// Difference across face
		scalar df = phiN - phiP;

		// Twice the differences across face-neighbour cells
		scalar tdcP = 2*(d & gradcP);
		scalar tdcN = 2*(d & gradcN);

		// Calculate the limiter according to the sign of the face difference
		scalar limiter;

		if (df > 0)
		{
			limiter = l_
				- k_*min(max(df - tdcP, 0), max(df - tdcN, 0))
				/(max(mag(df), max(mag(tdcP), mag(tdcN))) + SMALL);
		}
		else
		{
			limiter = l_
				- k_*min(max(tdcP - df, 0), max(tdcN - df, 0))
				/(max(mag(df), max(mag(tdcP), mag(tdcN))) + SMALL);
		}

		// Limit the limiter between linear and upwind
		return max(min(limiter, 1), 0);
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
