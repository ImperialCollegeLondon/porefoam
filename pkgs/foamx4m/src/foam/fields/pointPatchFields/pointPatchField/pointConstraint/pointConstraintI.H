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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pointConstraint::pointConstraint()
:
	Tuple2<label, vector>(0, vector::zero)
{}


inline Foam::pointConstraint::pointConstraint(Istream& is)
:
	Tuple2<label, vector>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointConstraint::applyConstraint(const vector& cd)
{
	if (first() == 0)
	{
		first() = 1;
		second() = cd;
	}
	else if (first() == 1)
	{
		vector planeNormal = cd ^ second();
		scalar magPlaneNormal = mag(planeNormal);

		if (magPlaneNormal > 1e-3)
		{
			first() = 2;
			second() = planeNormal/magPlaneNormal;
		}
	}
	else if (first() == 2)
	{
		if (mag(cd & second()) > 1e-3)
		{
			first() = 3;
			second() = vector::zero;
		}
	}
}


Foam::tensor Foam::pointConstraint::constraintTransformation() const
{
	if (first() == 0)
	{
		return I;
	}
	else if (first() == 1)
	{
		return I - sqr(second());
	}
	else if (first() == 2)
	{
		return sqr(second());
	}
	else
	{
		return tensor::zero;
	}
}


// ************************************************************************* //
