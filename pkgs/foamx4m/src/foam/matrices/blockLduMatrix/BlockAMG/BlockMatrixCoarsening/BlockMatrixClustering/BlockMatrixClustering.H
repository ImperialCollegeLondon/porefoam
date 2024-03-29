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
	BlockMatrixClustering

Description
	Block matrix AMG coarsening by Jasak clustering algorithm

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	BlockMatrixClustering.C

\*---------------------------------------------------------------------------*/

#ifndef BlockMatrixClustering_H
#define BlockMatrixClustering_H

#include "BlockMatrixCoarsening.H"
#include "BlockLduMatrix.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeff.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
class BlockMatrixClustering
:
	public BlockMatrixCoarsening<Type>
{
	// Private Data

		//- Reference to matrix
		const BlockLduMatrix<Type>& matrix_;

		//- Min group size
		const label minGroupSize_;

		//- Max group size
		const label maxGroupSize_;

		//- Reference to a templated norm calculator
		autoPtr<BlockCoeffNorm<Type> > normPtr_;

		//- Child array: for each fine equation give a clustering index
		labelList agglomIndex_;

		//- Face-restriction addressing
		mutable labelList coeffRestrictAddr_;

		//- Group size
		label groupSize_;

		//- Number of solo cells
		label nSolo_;

		//- Number of coarse equations
		label nCoarseEqns_;

		//- Can a coarse level be constructed?
		bool coarsen_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		BlockMatrixClustering(const BlockMatrixClustering<Type>&);

		// Disallow default bitwise assignment
		void operator=(const BlockMatrixClustering<Type>&);


		//- Calculate clustering index (child)
		void calcClustering();

		//- Restrict CoeffField.  Used for diag coefficient
		void restrictDiag
		(
			const CoeffField<Type>& Coeff,
			CoeffField<Type>& coarseCoeff
		) const;

		//- Agglomerate coeffs, symmetric matrix
		template<class DiagType, class ULType>
		void agglomerateCoeffs
		(
			Field<DiagType>& activeCoarseDiag,
			Field<ULType>& activeCoarseUpper,
			const Field<ULType>& activeFineUpper,
			const Field<ULType>& activeFineUpperTranspose
		) const;

		//- Agglomerate coeffs, assymmetric matrix
		template<class DiagType, class ULType>
		void agglomerateCoeffs
		(
			Field<DiagType>& activeCoarseDiag,
			Field<ULType>& activeCoarseUpper,
			const Field<ULType>& activeFineUpper,
			Field<ULType>& activeCoarseLower,
			const Field<ULType>& activeFineLower
		) const;

		//- Restrict CoeffField, decoupled version.  Used for diag coefficient
		void restrictDiagDecoupled
		(
			const CoeffField<Type>& Coeff,
			CoeffField<Type>& coarseCoeff
		) const;

	// Private Static Data

		//- Weighting factor
		static const debug::tolerancesSwitch weightFactor_;

		//- Diagonal scaling factor
		static const debug::tolerancesSwitch diagFactor_;


public:

	//- Runtime type information
	TypeName("AAMG");


	// Constructors

		//- Construct from matrix and group size
		BlockMatrixClustering
		(
			const BlockLduMatrix<Type>& matrix,
			const dictionary& dict,
			const label groupSize,
			const label minCoarseEqns
		);


	//- Destructor
	virtual ~BlockMatrixClustering();


	// Member Functions

		//- Can a coarse level be constructed?
		virtual bool coarsen() const
		{
			return coarsen_;
		}

		//- Restrict matrix
		virtual autoPtr<BlockAMGLevel<Type> > restrictMatrix() const;

		//- Restrict residual
		virtual void restrictResidual
		(
			const Field<Type>& res,
			Field<Type>& coarseRes
		) const;

		//- Prolongate correction
		virtual void prolongateCorrection
		(
			Field<Type>& x,
			const Field<Type>& coarseX
		) const;

		//- Update coarse matrix using same coefficients
		virtual void updateMatrix(BlockLduMatrix<Type>& coarseMatrix) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "BlockMatrixClustering.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
