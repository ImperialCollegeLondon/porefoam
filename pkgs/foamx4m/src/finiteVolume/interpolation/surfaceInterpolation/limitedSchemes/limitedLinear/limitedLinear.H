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
	Foam::limitedLinearLimiter

Description
	Class with limiter function which returns the limiter for the
	TVD limited linear differencing scheme based on r obtained from the
	LimiterFunc class.

	Used in conjunction with the template class LimitedScheme.

SourceFiles
	limitedLinear.C

\*---------------------------------------------------------------------------*/

#ifndef limitedLinear_H
#define limitedLinear_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class LimiterFunc>
class limitedLinearLimiter
:
	public LimiterFunc
{
	scalar k_;
	scalar twoByk_;

public:

	limitedLinearLimiter(Istream& is)
	:
		k_(readScalar(is))
	{
		if (k_ < 0 || k_ > 1)
		{
			FatalIOErrorIn("limitedLinearLimiter(Istream& is)", is)
				<< "coefficient = " << k_
				<< " should be >= 0 and <= 1"
				<< exit(FatalIOError);
		}

		// Avoid the /0 when k_ = 0
		twoByk_ = 2.0/max(k_, SMALL);
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
		scalar r = LimiterFunc::r
		(
			faceFlux, phiP, phiN, gradcP, gradcN, d
		);

		return max(min(twoByk_*r, 1), 0);
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
