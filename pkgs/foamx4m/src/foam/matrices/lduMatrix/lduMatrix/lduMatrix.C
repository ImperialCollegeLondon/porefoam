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

#include "lduMatrix.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(lduMatrix, 1);
}

const Foam::scalar Foam::lduMatrix::great_ = 1.0e+20;
const Foam::scalar Foam::lduMatrix::small_ = 1.0e-20;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::lduMatrix::lduMatrix(const lduMesh& mesh)
:
	lduMesh_(mesh),
	lowerPtr_(nullptr),
	diagPtr_(nullptr),
	upperPtr_(nullptr)
{}


Foam::lduMatrix::lduMatrix(const lduMatrix& A)
:
	lduMesh_(A.lduMesh_),
	lowerPtr_(nullptr),
	diagPtr_(nullptr),
	upperPtr_(nullptr)
{
	if (A.lowerPtr_)
	{
		lowerPtr_ = new scalarField(*(A.lowerPtr_));
	}

	if (A.diagPtr_)
	{
		diagPtr_ = new scalarField(*(A.diagPtr_));
	}

	if (A.upperPtr_)
	{
		upperPtr_ = new scalarField(*(A.upperPtr_));
	}
}


Foam::lduMatrix::lduMatrix(lduMatrix& A, bool reUse)
:
	lduMesh_(A.lduMesh_),
	lowerPtr_(nullptr),
	diagPtr_(nullptr),
	upperPtr_(nullptr)
{
	if (reUse)
	{
		if (A.lowerPtr_)
		{
			lowerPtr_ = A.lowerPtr_;
			A.lowerPtr_ = nullptr;
		}

		if (A.diagPtr_)
		{
			diagPtr_ = A.diagPtr_;
			A.diagPtr_ = nullptr;
		}

		if (A.upperPtr_)
		{
			upperPtr_ = A.upperPtr_;
			A.upperPtr_ = nullptr;
		}
	}
	else
	{
		if (A.lowerPtr_)
		{
			lowerPtr_ = new scalarField(*(A.lowerPtr_));
		}

		if (A.diagPtr_)
		{
			diagPtr_ = new scalarField(*(A.diagPtr_));
		}

		if (A.upperPtr_)
		{
			upperPtr_ = new scalarField(*(A.upperPtr_));
		}
	}
}


Foam::lduMatrix::lduMatrix
(
	const lduMesh& mesh,
	Istream& is
)
:
	lduMesh_(mesh),
	lowerPtr_(new scalarField(is)),
	diagPtr_(new scalarField(is)),
	upperPtr_(new scalarField(is))
{}


Foam::lduMatrix::~lduMatrix()
{
	if (lowerPtr_)
	{
		delete lowerPtr_;
	}

	if (diagPtr_)
	{
		delete diagPtr_;
	}

	if (upperPtr_)
	{
		delete upperPtr_;
	}
}


Foam::scalarField& Foam::lduMatrix::lower()
{
	if (!lowerPtr_)
	{
		if (upperPtr_)
		{
			lowerPtr_ = new scalarField(*upperPtr_);
		}
		else
		{
			lowerPtr_ = new scalarField(lduAddr().lowerAddr().size(), 0.0);
		}
	}

	return *lowerPtr_;
}


Foam::scalarField& Foam::lduMatrix::diag()
{
	if (!diagPtr_)
	{
		diagPtr_ = new scalarField(lduAddr().size(), 0.0);
	}

	return *diagPtr_;
}


Foam::scalarField& Foam::lduMatrix::upper()
{
	if (!upperPtr_)
	{
		if (lowerPtr_)
		{
			upperPtr_ = new scalarField(*lowerPtr_);
		}
		else
		{
			upperPtr_ = new scalarField(lduAddr().lowerAddr().size(), 0.0);
		}
	}

	return *upperPtr_;
}


const Foam::scalarField& Foam::lduMatrix::lower() const
{
	if (!lowerPtr_ && !upperPtr_)
	{
		FatalErrorIn("lduMatrix::lower() const")
			<< "lowerPtr_ or upperPtr_ unallocated"
			<< abort(FatalError);
	}

	if (lowerPtr_)
	{
		return *lowerPtr_;
	}
	else
	{
		return *upperPtr_;
	}
}


const Foam::scalarField& Foam::lduMatrix::diag() const
{
	if (!diagPtr_)
	{
		FatalErrorIn("const scalarField& lduMatrix::diag() const")
			<< "diagPtr_ unallocated"
			<< abort(FatalError);
	}

	return *diagPtr_;
}


const Foam::scalarField& Foam::lduMatrix::upper() const
{
	if (!lowerPtr_ && !upperPtr_)
	{
		FatalErrorIn("lduMatrix::upper() const")
			<< "lowerPtr_ or upperPtr_ unallocated"
			<< abort(FatalError);
	}

	if (upperPtr_)
	{
		return *upperPtr_;
	}
	else
	{
		return *lowerPtr_;
	}
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const lduMatrix& ldum)
{
	if (ldum.lowerPtr_)
	{
		os  << "Lower triangle = "
			<< *ldum.lowerPtr_
			<< endl << endl;
	}

	if (ldum.diagPtr_)
	{
		os  << "diagonal = "
			<< *ldum.diagPtr_
			<< endl << endl;
	}

	if (ldum.upperPtr_)
	{
		os  << "Upper triangle = "
			<< *ldum.upperPtr_
			<< endl << endl;
	}

	os.check("Ostream& operator<<(Ostream&, const lduMatrix&");

	return os;
}


// ************************************************************************* //
