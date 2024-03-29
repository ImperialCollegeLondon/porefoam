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

#include "procLduMatrix.H"
#include "procLduInterface.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::procLduMatrix::procLduMatrix
(
	const lduMatrix& ldum,
	const FieldField<Field, scalar>& interfaceCoeffs,
	const lduInterfaceFieldPtrsList& interfaces
)
:
	upperAddr_(ldum.lduAddr().upperAddr()),
	lowerAddr_(ldum.lduAddr().lowerAddr()),
	diag_(ldum.diag()),
	upper_(ldum.upper()),
	lower_(ldum.lower())
{
	label nInterfaces = 0;

	forAll(interfaces, i)
	{
		if (interfaces.set(i))
		{
			nInterfaces++;
		}
	}

	interfaces_.setSize(nInterfaces);

	nInterfaces = 0;

	forAll(interfaces, i)
	{
		if (interfaces.set(i))
		{
			interfaces_.set
			(
				nInterfaces++,
				new procLduInterface
				(
					interfaces[i],
					interfaceCoeffs[i]
				)
			);
		}
	}

}


Foam::procLduMatrix::procLduMatrix(Istream& is)
:
	upperAddr_(is),
	lowerAddr_(is),
	diag_(is),
	upper_(is),
	lower_(is),
	interfaces_(is)
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const procLduMatrix& cldum)
{
	os  << cldum.upperAddr_
		<< cldum.lowerAddr_
		<< cldum.diag_
		<< cldum.upper_
		<< cldum.lower_
		<< cldum.interfaces_;

	return os;
}


// ************************************************************************* //
