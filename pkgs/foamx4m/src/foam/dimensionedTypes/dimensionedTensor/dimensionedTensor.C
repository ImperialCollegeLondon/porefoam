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

#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
dimensionedTensor dimensionedTensor::T() const
{
	return dimensionedTensor
	(
		name()+".T()",
		dimensions(),
		value().T()
	);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

dimensionedScalar tr(const dimensionedTensor& dt)
{
	return dimensionedScalar
	(
		"tr("+dt.name()+')',
		dt.dimensions(),
		tr(dt.value())
	);
}


dimensionedTensor dev(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"dev("+dt.name()+')',
		dt.dimensions(),
		dev(dt.value())
	);
}


dimensionedTensor dev2(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"dev2("+dt.name()+')',
		dt.dimensions(),
		dev2(dt.value())
	);
}


dimensionedScalar det(const dimensionedTensor& dt)
{
	return dimensionedScalar
	(
		"det("+dt.name()+')',
		pow(dt.dimensions(), tensor::dim),
		det(dt.value())
	);
}


dimensionedTensor cof(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"cof("+dt.name()+')',
		dt.dimensions(),
		cof(dt.value())
	);
}


dimensionedTensor inv(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"inv("+dt.name()+')',
		dimless/dt.dimensions(),
		inv(dt.value())
	);
}


dimensionedTensor hinv(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"hinv("+dt.name()+')',
		dimless/dt.dimensions(),
		hinv(dt.value())
	);
}


dimensionedSymmTensor symm(const dimensionedTensor& dt)
{
	return dimensionedSymmTensor
	(
		"symm("+dt.name()+')',
		dt.dimensions(),
		symm(dt.value())
	);
}

dimensionedSymmTensor twoSymm(const dimensionedTensor& dt)
{
	return dimensionedSymmTensor
	(
		"twoSymm("+dt.name()+')',
		dt.dimensions(),
		twoSymm(dt.value())
	);
}

dimensionedTensor skew(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"skew("+dt.name()+')',
		dt.dimensions(),
		skew(dt.value())
	);
}


dimensionedVector eigenValues(const dimensionedTensor& dt)
{
	return dimensionedVector
	(
		"eigenValues("+dt.name()+')',
		dt.dimensions(),
		eigenValues(dt.value())
	);
}


dimensionedTensor eigenVectors(const dimensionedTensor& dt)
{
	return dimensionedTensor
	(
		"eigenVectors("+dt.name()+')',
		dimless,
		eigenVectors(dt.value())
	);
}


dimensionedVector eigenValues(const dimensionedSymmTensor& dt)
{
	return dimensionedVector
	(
		"eigenValues("+dt.name()+')',
		dt.dimensions(),
		eigenValues(dt.value())
	);
}


dimensionedTensor eigenVectors(const dimensionedSymmTensor& dt)
{
	return dimensionedTensor
	(
		"eigenVectors("+dt.name()+')',
		dimless,
		eigenVectors(dt.value())
	);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

dimensionedVector operator*(const dimensionedTensor& dt)
{
	return dimensionedVector
	(
		"*"+dt.name(),
		dt.dimensions(),
		*dt.value()
	);
}


dimensionedTensor operator*(const dimensionedVector& dv)
{
	return dimensionedTensor
	(
		"*"+dv.name(),
		dv.dimensions(),
		*dv.value()
	);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
