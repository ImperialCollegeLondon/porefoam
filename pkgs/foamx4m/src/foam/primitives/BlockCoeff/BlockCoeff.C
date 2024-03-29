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

Description
	Generic block coefficient type.  Used in BlockLduMatrix.  HJ, 2/Jan/2006

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
typename Foam::BlockCoeff<Type>::scalarType&
Foam::BlockCoeff<Type>::toScalar()
{
	if (!scalarCoeffPtr_)
	{
		// Debug check: demotion
		if (linearCoeffPtr_ || squareCoeffPtr_)
		{
			FatalErrorIn
			(
				"BlockCoeff<Type>::scalarType& "
				"BlockCoeff<Type>::toScalar()"
			)   << "Detected demotion to scalar.  Probably an error"
				<< abort(FatalError);
		}

		scalarCoeffPtr_ = new scalarType(pTraits<scalarType>::zero);
	}

	return *scalarCoeffPtr_;
}


template<class Type>
typename Foam::BlockCoeff<Type>::linearType&
Foam::BlockCoeff<Type>::toLinear()
{
	if (!linearCoeffPtr_)
	{
		// Debug check: demotion
		if (squareCoeffPtr_)
		{
			FatalErrorIn
			(
				"BlockCoeff<Type>::linearType& "
				"BlockCoeff<Type>::toLinear()"
			)   << "Detected demotion to linear.  Probably an error"
				<< abort(FatalError);
		}

		linearCoeffPtr_ = new linearType(pTraits<linearType>::zero);

		// If scalar is active, promote to linear
		if (scalarCoeffPtr_)
		{
			*linearCoeffPtr_ = (*scalarCoeffPtr_)*pTraits<linearType>::one;
			deleteDemandDrivenData(scalarCoeffPtr_);
		}
	}

	return *linearCoeffPtr_;
}


template<class Type>
typename Foam::BlockCoeff<Type>::squareType&
Foam::BlockCoeff<Type>::toSquare()
{
	if (!squareCoeffPtr_)
	{
		squareCoeffPtr_ = new squareType(pTraits<squareType>::zero);

		// If scalar is active, promote to square
		if (scalarCoeffPtr_)
		{
			expandScalar(*squareCoeffPtr_, *scalarCoeffPtr_);
			deleteDemandDrivenData(scalarCoeffPtr_);
		}

		// If linear is active, promote to square
		if (linearCoeffPtr_)
		{
			expandLinear(*squareCoeffPtr_, *linearCoeffPtr_);
			deleteDemandDrivenData(linearCoeffPtr_);
		}
	}

	return *squareCoeffPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeff<Type>::BlockCoeff()
:
	scalarCoeffPtr_(nullptr),
	linearCoeffPtr_(nullptr),
	squareCoeffPtr_(nullptr)
{}


template<class Type>
Foam::BlockCoeff<Type>::BlockCoeff(const BlockCoeff<Type>& f)
:
	scalarCoeffPtr_(nullptr),
	linearCoeffPtr_(nullptr),
	squareCoeffPtr_(nullptr)
{
	if (f.scalarCoeffPtr_)
	{
		scalarCoeffPtr_ = new scalarType(*(f.scalarCoeffPtr_));
	}
	else if (f.linearCoeffPtr_)
	{
		linearCoeffPtr_ = new linearType(*(f.linearCoeffPtr_));
	}
	else if (f.squareCoeffPtr_)
	{
		squareCoeffPtr_ = new squareType(*(f.squareCoeffPtr_));
	}
}


template<class Type>
Foam::BlockCoeff<Type>::BlockCoeff(Istream& is)
:
	scalarCoeffPtr_(nullptr),
	linearCoeffPtr_(nullptr),
	squareCoeffPtr_(nullptr)
{
	// Read keyword and pick up allocated field
	word key(is);

	if
	(
		key
	 == blockCoeffBase::activeLevelNames_[blockCoeffBase::UNALLOCATED]
	)
	{
	}
	else if
	(
		key
	 == blockCoeffBase::activeLevelNames_[blockCoeffBase::SCALAR]
	)
	{
		scalarCoeffPtr_ = new scalarType(readScalar(is));
	}
	else if
	(
		key
	 == blockCoeffBase::activeLevelNames_[blockCoeffBase::LINEAR]
	)
	{
		linearCoeffPtr_ = new linearType(is);
	}
	else if
	(
		key
	 == blockCoeffBase::activeLevelNames_[blockCoeffBase::SQUARE]
	)
	{
		squareCoeffPtr_ = new squareType(is);
	}
	else
	{
		FatalIOErrorIn
		(
			"BlockCoeff<Type>::BlockCoeff(Istream& is)",
			is
		)   << "invalid keyword while reading: " << key
			<< exit(FatalIOError);
	}
}


template<class Type>
Foam::BlockCoeff<Type> Foam::BlockCoeff<Type>::clone() const
{
	return BlockCoeff<Type>(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockCoeff<Type>::~BlockCoeff()
{
	this->clear();
}


template<class Type>
inline void Foam::BlockCoeff<Type>::clear()
{
	deleteDemandDrivenData(scalarCoeffPtr_);
	deleteDemandDrivenData(linearCoeffPtr_);
	deleteDemandDrivenData(squareCoeffPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::blockCoeffBase::activeLevel
Foam::BlockCoeff<Type>::activeType() const
{
	if (scalarCoeffPtr_)
	{
		return blockCoeffBase::SCALAR;
	}
	else if (linearCoeffPtr_)
	{
		return blockCoeffBase::LINEAR;
	}
	else if (squareCoeffPtr_)
	{
		return blockCoeffBase::SQUARE;
	}
	else
	{
		return blockCoeffBase::UNALLOCATED;
	}
}


template<class Type>
void Foam::BlockCoeff<Type>::checkActive() const
{
	label nActive = 0;

	if (scalarCoeffPtr_) nActive++;
	if (linearCoeffPtr_) nActive++;
	if (squareCoeffPtr_) nActive++;

	if (nActive > 1)
	{
		FatalErrorIn("void Foam::BlockCoeff<Type>::checkActive() const")
			<< "Activation/deactivation error.  nActive = " << nActive
			<< abort(FatalError);
	}
}


template<class Type>
const typename Foam::BlockCoeff<Type>::scalarType&
Foam::BlockCoeff<Type>::asScalar() const
{
	if (!scalarCoeffPtr_)
	{
		FatalErrorIn
		(
			"BlockCoeff<Type>::scalarType& "
			"BlockCoeff<Type>::asScalar()"
		)   << "Requested scalar but active type is: "
			<< blockCoeffBase::activeLevelNames_[this->activeType()]
			<< ".  This is not allowed."
			<< abort(FatalError);
	}

	return *scalarCoeffPtr_;
}


template<class Type>
const typename Foam::BlockCoeff<Type>::linearType&
Foam::BlockCoeff<Type>::asLinear() const
{
	if (!linearCoeffPtr_)
	{
		FatalErrorIn
		(
			"BlockCoeff<Type>::linearType& "
			"BlockCoeff<Type>::asLinear()"
		)   << "Requested linear but active type is: "
			<< blockCoeffBase::activeLevelNames_[this->activeType()]
			<< ".  This is not allowed."
			<< abort(FatalError);
	}

	return *linearCoeffPtr_;
}


template<class Type>
const typename Foam::BlockCoeff<Type>::squareType&
Foam::BlockCoeff<Type>::asSquare() const
{
	if (!squareCoeffPtr_)
	{
		FatalErrorIn
		(
			"BlockCoeff<Type>::squareType& "
			"BlockCoeff<Type>::asSquare()"
		)   << "Requested square but active type is: "
			<< blockCoeffBase::activeLevelNames_[this->activeType()]
			<< ".  This is not allowed."
			<< abort(FatalError);
	}

	return *squareCoeffPtr_;
}


template<class Type>
typename Foam::BlockCoeff<Type>::scalarType&
Foam::BlockCoeff<Type>::asScalar()
{
	if (linearCoeffPtr_ || squareCoeffPtr_)
	{
		FatalErrorIn
		(
			"BlockCoeff<Type>::scalarType& "
			"BlockCoeff<Type>::asScalar()"
		)   << "Requested scalar but active type is: "
			<< blockCoeffBase::activeLevelNames_[this->activeType()]
			<< ".  This is not allowed."
			<< abort(FatalError);
	}

	if (!scalarCoeffPtr_)
	{
		return this->toScalar();
	}

	return *scalarCoeffPtr_;
}


template<class Type>
typename Foam::BlockCoeff<Type>::linearType&
Foam::BlockCoeff<Type>::asLinear()
{
	if (squareCoeffPtr_)
	{
		FatalErrorIn
		(
			"BlockCoeff<Type>::linearType& "
			"BlockCoeff<Type>::asLinear()"
		)   << "Requested linear but active type is: "
			<< blockCoeffBase::activeLevelNames_[this->activeType()]
			<< ".  This is not allowed."
			<< abort(FatalError);
	}

	if (!linearCoeffPtr_)
	{
		return this->toLinear();
	}

	return *linearCoeffPtr_;
}


template<class Type>
typename Foam::BlockCoeff<Type>::squareType&
Foam::BlockCoeff<Type>::asSquare()
{
	if (!squareCoeffPtr_)
	{
		return this->toSquare();
	}

	return *squareCoeffPtr_;
}


template<class Type>
typename Foam::BlockCoeff<Type>::scalarType
Foam::BlockCoeff<Type>::component(const direction dir) const
{
	if (scalarCoeffPtr_)
	{
		return *scalarCoeffPtr_;
	}
	else if (linearCoeffPtr_)
	{
		return linearCoeffPtr_->component(dir);
	}
	else if (squareCoeffPtr_)
	{
		return contractLinear
		(
			*squareCoeffPtr_
		).component(dir);
	}
	else
	{
		FatalErrorIn
		(
			"tmp<BlockCoeff<Type>::scalarType>"
			"BlockCoeff<Type>::component(const direction dir) const"
		)   << " not allocated."
			<< abort(FatalError);
	}

	// Dummy return to keep compiler happy
	return *scalarCoeffPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockCoeff<Type>::operator=(const BlockCoeff<Type>& f)
{
	if (this == &f)
	{
		FatalErrorIn("BlockCoeff<Type>::operator=(const BlockCoeff<Type>&)")
			<< "attempted assignment to self"
			<< abort(FatalError);
	}

	if (f.scalarCoeffPtr_)
	{
		this->toScalar() = *(f.scalarCoeffPtr_);
	}
	else if (f.linearCoeffPtr_)
	{
		this->toLinear() = *(f.linearCoeffPtr_);
	}
	else if (f.squareCoeffPtr_)
	{
		this->toSquare() = *(f.squareCoeffPtr_);
	}
	else
	{
		// Not allocated - do nothing
	}
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const BlockCoeff<Type>& f)
{
	// Write active type
	os << blockCoeffBase::activeLevelNames_[f.activeType()] << nl;

	if (f.scalarCoeffPtr_)
	{
		os << *(f.scalarCoeffPtr_);
	}
	else if (f.linearCoeffPtr_)
	{
		os << *(f.linearCoeffPtr_);
	}
	else if (f.squareCoeffPtr_)
	{
		os << *(f.squareCoeffPtr_);
	}

	return os;
}


// ************************************************************************* //
