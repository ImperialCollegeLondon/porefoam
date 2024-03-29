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

#include "tensor2D.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const tensor2D::typeName = "tensor2D";

template<>
const char* tensor2D::componentNames[] =
{
	"xx", "xy",
	"yx", "yy"
};

template<>
const tensor2D tensor2D::zero
(
	0, 0,
	0, 0
);

template<>
const tensor2D tensor2D::one
(
	1, 1,
	1, 1
);

template<>
const tensor2D tensor2D::max
(
	VGREAT, VGREAT,
	VGREAT, VGREAT
);

template<>
const tensor2D tensor2D::min
(
	-VGREAT, -VGREAT,
	-VGREAT, -VGREAT
);

template<>
const tensor2D tensor2D::I
(
	1, 0,
	0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return eigenvalues in ascending order of absolute values
vector2D eigenValues(const tensor2D& t)
{
	scalar i = 0;
	scalar ii = 0;

	if (mag(t.xy()) < SMALL && mag(t.yx()) < SMALL)
	{
		i = t.xx();
		ii = t.yy();
	}
	else
	{
		scalar mb = t.xx() + t.yy();
		scalar c = t.xx()*t.yy() - t.xy()*t.yx();

		// If there is a zero root
		if (mag(c) < SMALL)
		{
			i = 0;
			ii = mb;
		}
		else
		{
			scalar disc = sqr(mb) - 4*c;

			if (disc > 0)
			{
				scalar q = sqrt(disc);

				i = 0.5*(mb - q);
				ii = 0.5*(mb + q);
			}
			else
			{
				FatalErrorIn("eigenValues(const tensor2D&)")
					<< "zero and complex eigenvalues in tensor2D: " << t
					<< abort(FatalError);
			}
		}
	}

	// Sort the eigenvalues into ascending order
	if (mag(i) > mag(ii))
	{
		Swap(i, ii);
	}

	return vector2D(i, ii);
}


vector2D eigenVector(const tensor2D& t, const scalar lambda)
{
	if (lambda < SMALL)
	{
		return vector2D::zero;
	}

	if (mag(t.xy()) < SMALL && mag(t.yx()) < SMALL)
	{
		if (lambda > min(t.xx(), t.yy()))
		{
			return vector2D(1, 0);
		}
		else
		{
			return vector2D(0, 1);
		}
	}
	else if (mag(t.xy()) < SMALL)
	{
		return vector2D(lambda - t.yy(), t.yx());
	}
	else
	{
		return vector2D(t.xy(), lambda - t.yy());
	}
}


tensor2D eigenVectors(const tensor2D& t)
{
	vector2D evals(eigenValues(t));

	tensor2D evs
	(
		eigenVector(t, evals.x()),
		eigenVector(t, evals.y())
	);

	return evs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
