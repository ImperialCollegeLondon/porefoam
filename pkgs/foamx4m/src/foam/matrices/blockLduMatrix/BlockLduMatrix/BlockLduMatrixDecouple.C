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
	Vector-matrix multiplication operations for a block matrix

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledSumDiag()
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	TypeCoeffField& Diag = this->diag();

	const unallocLabelList& l = lduAddr().lowerAddr();
	const unallocLabelList& u = lduAddr().upperAddr();

	if (this->symmetric())
	{
		// Symmetric matrix: re-use upper for lower coefficients

		const TypeCoeffField& Upper =
			const_cast<const BlockLduMatrix<Type>&>(*this).upper();

		if
		(
			Upper.activeType() == blockCoeffBase::LINEAR
		 || Diag.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiag = Diag.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] += activeUpper[coeffI];
				activeDiag[u[coeffI]] += activeUpper[coeffI];
			}
		}
		else if
		(
			Upper.activeType() == blockCoeffBase::SCALAR
		 || Diag.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiag = Diag.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] += activeUpper[coeffI];
				activeDiag[u[coeffI]] += activeUpper[coeffI];
			}
		}
	}
	else if (this->asymmetric())
	{
		// Full asymmetric matrix

		const TypeCoeffField& Lower =
			const_cast<const BlockLduMatrix<Type>&>(*this).lower();

		const TypeCoeffField& Upper =
			const_cast<const BlockLduMatrix<Type>&>(*this).upper();

		if
		(
			Lower.activeType() == blockCoeffBase::LINEAR
		 || Upper.activeType() == blockCoeffBase::LINEAR
		 || Diag.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeLower = Lower.asLinear();
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiag = Diag.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] += activeLower[coeffI];
				activeDiag[u[coeffI]] += activeUpper[coeffI];
			}
		}
		else if
		(
			Lower.activeType() == blockCoeffBase::SCALAR
		 || Upper.activeType() == blockCoeffBase::SCALAR
		 || Diag.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeLower = Lower.asScalar();
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiag = Diag.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] += activeLower[coeffI];
				activeDiag[u[coeffI]] += activeUpper[coeffI];
			}
		}
	}
	else
	{
		FatalErrorIn("void BlockLduMatrix<Type>::decoupledSumDiag()")
			<< "No off-diagonal available"
			<< abort(FatalError);
	}
}


template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledNegSumDiag()
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	TypeCoeffField& Diag = this->diag();

	const unallocLabelList& l = lduAddr().lowerAddr();
	const unallocLabelList& u = lduAddr().upperAddr();

	if (this->symmetric())
	{
		// Symmetric matrix: re-use upper for lower coefficients

		const TypeCoeffField& Upper =
			const_cast<const BlockLduMatrix<Type>&>(*this).upper();

		if
		(
			Upper.activeType() == blockCoeffBase::LINEAR
		 || Diag.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiag = Diag.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] -= activeUpper[coeffI];
				activeDiag[u[coeffI]] -= activeUpper[coeffI];
			}
		}
		else if
		(
			Upper.activeType() == blockCoeffBase::SCALAR
		 || Diag.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiag = Diag.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] -= activeUpper[coeffI];
				activeDiag[u[coeffI]] -= activeUpper[coeffI];
			}
		}
	}
	else if (this->asymmetric())
	{
		// Full asymmetric matrix

		const TypeCoeffField& Lower =
			const_cast<const BlockLduMatrix<Type>&>(*this).lower();

		const TypeCoeffField& Upper =
			const_cast<const BlockLduMatrix<Type>&>(*this).upper();

		if
		(
			Lower.activeType() == blockCoeffBase::LINEAR
		 || Upper.activeType() == blockCoeffBase::LINEAR
		 || Diag.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeLower = Lower.asLinear();
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiag = Diag.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] -= activeLower[coeffI];
				activeDiag[u[coeffI]] -= activeUpper[coeffI];
			}
		}
		else if
		(
			Lower.activeType() == blockCoeffBase::SCALAR
		 || Upper.activeType() == blockCoeffBase::SCALAR
		 || Diag.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeLower = Lower.asScalar();
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiag = Diag.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiag[l[coeffI]] -= activeLower[coeffI];
				activeDiag[u[coeffI]] -= activeUpper[coeffI];
			}
		}
	}
	else
	{
		FatalErrorIn("void BlockLduMatrix<Type>::decoupledNegSumDiag()")
			<< "No off-diagonal available"
			<< abort(FatalError);
	}
}


template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledCheck() const
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	// Copy the diagonal
	TypeCoeffField DiagCopy(this->diag().size());

	const unallocLabelList& l = lduAddr().lowerAddr();
	const unallocLabelList& u = lduAddr().upperAddr();

	if (this->symmetric())
	{
		// Symmetric matrix: re-use upper for lower coefficients

		const TypeCoeffField& Upper = this->upper();

		if
		(
			Upper.activeType() == blockCoeffBase::LINEAR
		 || DiagCopy.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiagCopy = DiagCopy.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiagCopy[l[coeffI]] += activeUpper[coeffI];
				activeDiagCopy[u[coeffI]] += activeUpper[coeffI];
			}

			Info<< "void BlockLduMatrix<Type>::decoupledCheck() const : "
				<< "Symmetric matrix: raw matrix difference: "
				<< sum(mag(activeDiagCopy))
				<< " scaled: "
				<< sum(mag(activeDiagCopy))/sum(mag(this->diag().asLinear()))
				<< endl;
		}
		else if
		(
			Upper.activeType() == blockCoeffBase::SCALAR
		 || DiagCopy.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiagCopy = DiagCopy.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiagCopy[l[coeffI]] += activeUpper[coeffI];
				activeDiagCopy[u[coeffI]] += activeUpper[coeffI];
			}

			Info<< "void BlockLduMatrix<Type>::decoupledCheck() const : "
				<< "Symmetric matrix: raw matrix difference: "
				<< sum(mag(activeDiagCopy))
				<< " scaled: "
				<< sum(mag(activeDiagCopy))/sum(mag(this->diag().asScalar()))
				<< endl;
		}
	}
	else if (this->asymmetric())
	{
		// Full asymmetric matrix

		const TypeCoeffField& Lower = this->lower();
		const TypeCoeffField& Upper = this->upper();

		if
		(
			Lower.activeType() == blockCoeffBase::LINEAR
		 || Upper.activeType() == blockCoeffBase::LINEAR
		 || DiagCopy.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeLower = Lower.asLinear();
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiagCopy = DiagCopy.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiagCopy[l[coeffI]] += activeLower[coeffI];
				activeDiagCopy[u[coeffI]] += activeUpper[coeffI];
			}

			Info<< "void BlockLduMatrix<Type>::decoupledCheck() const : "
				<< "Asymmetric matrix: raw matrix difference: "
				<< sum(mag(activeDiagCopy))
				<< " scaled: "
				<< sum(mag(activeDiagCopy))/sum(mag(this->diag().asLinear()))
				<< endl;
		}
		else if
		(
			Lower.activeType() == blockCoeffBase::SCALAR
		 || Upper.activeType() == blockCoeffBase::SCALAR
		 || DiagCopy.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeLower = Lower.asScalar();
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiagCopy = DiagCopy.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeDiagCopy[l[coeffI]] += activeLower[coeffI];
				activeDiagCopy[u[coeffI]] += activeUpper[coeffI];
			}

			Info<< "void BlockLduMatrix<Type>::decoupledCheck() const : "
				<< "Asymmetric matrix: raw matrix difference: "
				<< sum(mag(activeDiagCopy))
				<< " scaled: "
				<< sum(mag(activeDiagCopy))/sum(mag(this->diag().asScalar()))
				<< endl;
		}
	}
	else
	{
		Info<< "void BlockLduMatrix<Type>::decoupledCheck() const : "
			<< "Diagonal matrix" << endl;
	}
}


template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledRelax
(
	const TypeField& x,
	TypeField& b,
	const scalar alpha
)
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	//HJ Missing code: add coupling coefficients to under-relaxation
	//   HJ, 21/Feb/2008

	if (alpha <= 0)
	{
		return;
	}

	TypeCoeffField& Diag = this->diag();

	// Create multiplication function object
	typename BlockCoeff<Type>::multiply mult;

	const unallocLabelList& l = lduAddr().lowerAddr();
	const unallocLabelList& u = lduAddr().upperAddr();

	if (this->symmetric())
	{
		// Symmetric matrix: re-use upper for lower coefficients

		const TypeCoeffField& Upper =
			const_cast<const BlockLduMatrix<Type>&>(*this).upper();

		if
		(
			Upper.activeType() == blockCoeffBase::LINEAR
		 || Diag.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiag = Diag.asLinear();

			// Make a copy of diagonal before relaxation
			linearTypeField activeDiagOld = activeDiag;

			linearTypeField sumOff
			(
				activeDiag.size(),
				pTraits<typename TypeCoeffField::linearType>::zero
			);

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				sumOff[u[coeffI]] += cmptMag(activeUpper[coeffI]);
				sumOff[l[coeffI]] += cmptMag(activeUpper[coeffI]);
			}

			activeDiag = max(activeDiag, sumOff);
			activeDiag *= 1.0/alpha;

			// Add the relaxation contribution to b
			forAll (b, i)
			{
				b[i] += mult(activeDiag[i] - activeDiagOld[i], x[i]);
			}
		}
		else if
		(
			Upper.activeType() == blockCoeffBase::SCALAR
		 || Diag.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiag = Diag.asScalar();

			// Make a copy of diagonal before relaxation
			scalarTypeField activeDiagOld = activeDiag;

			scalarTypeField sumOff
			(
				activeDiag.size(),
				pTraits<typename TypeCoeffField::scalarType>::zero
			);

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				sumOff[u[coeffI]] += mag(activeUpper[coeffI]);
				sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
			}

			activeDiag = max(activeDiag, sumOff);
			activeDiag *= 1.0/alpha;

			// Add the relaxation contribution to b
			forAll (b, i)
			{
				b[i] += mult(activeDiag[i] - activeDiagOld[i], x[i]);
			}
		}
	}
	else if (this->asymmetric())
	{
		// Full asymmetric matrix

		const TypeCoeffField& Lower =
			const_cast<const BlockLduMatrix<Type>&>(*this).lower();

		const TypeCoeffField& Upper =
			const_cast<const BlockLduMatrix<Type>&>(*this).upper();

		if
		(
			Lower.activeType() == blockCoeffBase::LINEAR
		 || Upper.activeType() == blockCoeffBase::LINEAR
		 || Diag.activeType() == blockCoeffBase::LINEAR
		)
		{
			const linearTypeField& activeLower = Lower.asLinear();
			const linearTypeField& activeUpper = Upper.asLinear();
			linearTypeField& activeDiag = Diag.asLinear();

			linearTypeField sumOff
			(
				activeDiag.size(),
				pTraits<typename TypeCoeffField::linearType>::zero
			);

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				sumOff[u[coeffI]] += cmptMag(activeLower[coeffI]);
				sumOff[l[coeffI]] += cmptMag(activeUpper[coeffI]);
			}

			activeDiag = max(activeDiag, sumOff);
			activeDiag *= 1.0/alpha;
		}
		else if
		(
			Lower.activeType() == blockCoeffBase::SCALAR
		 || Upper.activeType() == blockCoeffBase::SCALAR
		 || Diag.activeType() == blockCoeffBase::SCALAR
		)
		{
			const scalarTypeField& activeLower = Lower.asScalar();
			const scalarTypeField& activeUpper = Upper.asScalar();
			scalarTypeField& activeDiag = Diag.asScalar();

			// Make a copy of diagonal before relaxation
			scalarTypeField activeDiagOld = activeDiag;

			scalarTypeField sumOff
			(
				activeDiag.size(),
				pTraits<typename TypeCoeffField::scalarType>::zero
			);

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				sumOff[u[coeffI]] += mag(activeLower[coeffI]);
				sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
			}

			activeDiag = max(activeDiag, sumOff);
			activeDiag *= 1.0/alpha;

			forAll (b, i)
			{
				b[i] += (activeDiag[i] - activeDiagOld[i])*x[i];
			}
		}
	}
}


template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledMultEqOp(const scalarField& sf)
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	if (diagPtr_)
	{
		*diagPtr_ *= sf;
	}

	if (upperPtr_)
	{
		TypeCoeffField& Upper = *upperPtr_;

		const unallocLabelList& l = lduAddr().lowerAddr();

		if (Upper.activeType() == blockCoeffBase::SCALAR)
		{
			scalarTypeField& activeUpper = Upper.asScalar();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeUpper[coeffI] *= sf[l[coeffI]];
			}
		}
		else if (Upper.activeType() == blockCoeffBase::LINEAR)
		{
			linearTypeField& activeUpper = Upper.asLinear();

			for (label coeffI = 0; coeffI < l.size(); coeffI++)
			{
				activeUpper[coeffI] *= sf[l[coeffI]];
			}
		}
	}

	if (lowerPtr_)
	{
		TypeCoeffField& Lower = *lowerPtr_;

		const unallocLabelList& u = lduAddr().upperAddr();

		if (Lower.activeType() == blockCoeffBase::SCALAR)
		{
			scalarTypeField& activeLower = Lower.asScalar();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				activeLower[coeffI] *= sf[u[coeffI]];
			}
		}
		else if (Lower.activeType() == blockCoeffBase::LINEAR)
		{
			linearTypeField& activeLower = Lower.asLinear();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				activeLower[coeffI] *= sf[u[coeffI]];
			}
		}
	}
}


template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledAmulCore
(
	TypeField& Ax,
	const TypeField& x
) const
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	const unallocLabelList& u = lduAddr().upperAddr();
	const unallocLabelList& l = lduAddr().lowerAddr();

	// In order to do automatic multiplication, diagonal needs to be recognised
	// as a decoupled coeff field.  HJ, 19/Feb/2008
//     const TypeCoeffField& Diag = this->diag();
	const DecoupledCoeffField<Type>& Diag = this->diag();
	const TypeCoeffField& Upper = this->upper();

	// Diagonal multiplication, no indirection
	multiply(Ax, Diag, x);

	// Create multiplication function object
	typename BlockCoeff<Type>::multiply mult;

	// Lower multiplication

	if (symmetric())
	{
		if (Upper.activeType() == blockCoeffBase::SCALAR)
		{
			const scalarTypeField& activeUpper = Upper.asScalar();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Ax[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
			}
		}
		else if (Upper.activeType() == blockCoeffBase::LINEAR)
		{
			const linearTypeField& activeUpper = Upper.asLinear();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Ax[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
			}
		}
	}
	else // Asymmetric matrix
	{
		const TypeCoeffField& Lower = this->lower();

		if (Lower.activeType() == blockCoeffBase::SCALAR)
		{
			const scalarTypeField& activeLower = Lower.asScalar();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Ax[u[coeffI]] += mult(activeLower[coeffI], x[l[coeffI]]);
			}
		}
		else if (Lower.activeType() == blockCoeffBase::LINEAR)
		{
			const linearTypeField& activeLower = Lower.asLinear();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Ax[u[coeffI]] += mult(activeLower[coeffI], x[l[coeffI]]);
			}
		}
	}


	// Upper multiplication

	if (Upper.activeType() == blockCoeffBase::SCALAR)
	{
		const scalarTypeField& activeUpper = Upper.asScalar();

		for (label coeffI = 0; coeffI < u.size(); coeffI++)
		{
			Ax[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
		}
	}
	else if (Upper.activeType() == blockCoeffBase::LINEAR)
	{
		const linearTypeField& activeUpper = Upper.asLinear();

		for (label coeffI = 0; coeffI < u.size(); coeffI++)
		{
			Ax[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
		}
	}
}


template<class Type>
void Foam::BlockLduMatrix<Type>::decoupledTmulCore
(
	TypeField& Tx,
	const TypeField& x
) const
{
	typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
	typedef typename TypeCoeffField::linearTypeField linearTypeField;

	const unallocLabelList& u = lduAddr().upperAddr();
	const unallocLabelList& l = lduAddr().lowerAddr();

	// In order to do automatic multiplication, diagonal needs to be recognised
	// as a decoupled coeff field.  HJ, 19/Feb/2008
	const DecoupledCoeffField<Type>& Diag = this->diag();
	const TypeCoeffField& Upper = this->upper();

	// Diagonal multiplication, no indirection
	multiply(Tx, Diag, x);

	// Create multiplication function object
	typename BlockCoeff<Type>::multiply mult;

	// Upper multiplication

	if (Upper.activeType() == blockCoeffBase::SCALAR)
	{
		const scalarTypeField& activeUpper = Upper.asScalar();

		for (label coeffI = 0; coeffI < u.size(); coeffI++)
		{
			Tx[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
		}
	}
	else if (Upper.activeType() == blockCoeffBase::LINEAR)
	{
		const linearTypeField& activeUpper = Upper.asLinear();

		for (label coeffI = 0; coeffI < u.size(); coeffI++)
		{
			Tx[u[coeffI]] += mult(activeUpper[coeffI], x[l[coeffI]]);
		}
	}

	// Lower multiplication

	if (symmetric())
	{
		if (Upper.activeType() == blockCoeffBase::SCALAR)
		{
			const scalarTypeField& activeUpper = Upper.asScalar();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Tx[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
			}
		}
		else if (Upper.activeType() == blockCoeffBase::LINEAR)
		{
			const linearTypeField& activeUpper = Upper.asLinear();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Tx[l[coeffI]] += mult(activeUpper[coeffI], x[u[coeffI]]);
			}
		}
	}
	else // Asymmetric matrix
	{
		const TypeCoeffField& Lower = this->lower();

		if (Lower.activeType() == blockCoeffBase::SCALAR)
		{
			const scalarTypeField& activeLower = Lower.asScalar();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Tx[l[coeffI]] += mult(activeLower[coeffI], x[u[coeffI]]);
			}
		}
		else if (Lower.activeType() == blockCoeffBase::LINEAR)
		{
			const linearTypeField& activeLower = Lower.asLinear();

			for (label coeffI = 0; coeffI < u.size(); coeffI++)
			{
				Tx[l[coeffI]] += mult(activeLower[coeffI], x[u[coeffI]]);
			}
		}
	}
}


// ************************************************************************* //
