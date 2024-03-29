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

#include "fvBlockMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvBlockMatrix<Type>::negate()
{
	BlockLduSystem<Type, Type>::negate();
	psi_.negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvBlockMatrix<Type>::operator=
(
	const fvBlockMatrix<Type>& bxs
)
{
	if (this == &bxs)
	{
		FatalErrorIn
		(
			"void fvBlockMatrix<Type>::operator="
			"(const fvBlockMatrix<Type>& bs)"
		)   << "attempted assignment to self"
			<< abort(FatalError);
	}

	BlockLduSystem<Type, Type>::operator=(bxs);
	psi_ = bxs.psi();
}


template<class Type>
void Foam::fvBlockMatrix<Type>::operator+=
(
	const fvBlockMatrix<Type>& bxs
)
{
	BlockLduSystem<Type, Type>::operator+=(bxs);
	psi_ += bxs.psi();
}

template<class Type>
void Foam::fvBlockMatrix<Type>::operator-=
(
	const fvBlockMatrix<Type>& bxs
)
{
	BlockLduSystem<Type, Type>::operator-=(bxs);
	psi_ -= bxs.psi();
}

template<class Type>
void Foam::fvBlockMatrix<Type>::operator*=
(
	const scalarField& sf
)
{
	BlockLduSystem<Type, Type>::operator*=(sf);
	psi_ *= sf;
}

template<class Type>
void Foam::fvBlockMatrix<Type>::operator*=
(
	const scalar s
)
{
	BlockLduSystem<Type, Type>::operator*=(s);
	psi_ *= s;
}


// ************************************************************************* //
