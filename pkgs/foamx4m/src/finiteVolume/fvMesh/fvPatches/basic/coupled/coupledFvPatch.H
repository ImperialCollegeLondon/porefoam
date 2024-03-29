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
	Foam::coupledFvPatch

Description
	An abstract base class for patches that couple regions of the
	computational domain e.g. cyclic and processor-processor links.

SourceFiles
	coupledFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFvPatch_H
#define coupledFvPatch_H

#include "fvPatch.H"
#include "lduInterface.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class coupledFvPatch
:
	public lduInterface,
	public fvPatch
{
	// Private data

		const coupledPolyPatch& coupledPolyPatch_;


protected:

	// Protected Member Functions

		//- Make patch weighting factors
		virtual void makeWeights(fvsPatchScalarField&) const = 0;

		//- Make patch face - neighbour cell distances
		virtual void makeDeltaCoeffs(fvsPatchScalarField&) const = 0;

		//- Make patch face non-orthogonality correction vectors
		virtual void makeCorrVecs(fvsPatchVectorField&) const;


public:

	friend class surfaceInterpolation;


	//- Runtime type information
	TypeName(coupledPolyPatch::typeName_());


	// Constructors

		//- Construct from polyPatch
		coupledFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
		:
			fvPatch(patch, bm),
			coupledPolyPatch_(refCast<const coupledPolyPatch>(patch))
		{}


	//- Destructor
	virtual ~coupledFvPatch();


	// Member Functions

		// Access

			//- Return true if patch is coupled
			virtual bool coupled() const
			{
				return coupledPolyPatch_.coupled();
			}

			//- Return true if patch is coupled
			virtual bool master() const
			{
				return coupledPolyPatch_.master();
			}

			//- Return face transformation tensor
			const tensorField& forwardT() const
			{
				return coupledPolyPatch_.forwardT();
			}

			//- Return neighbour-cell transformation tensor
			const tensorField& reverseT() const
			{
				return coupledPolyPatch_.reverseT();
			}

			//- Are the cyclic planes parallel
			bool parallel() const
			{
				return coupledPolyPatch_.parallel();
			}

			//- Return faceCell addressing
			virtual const unallocLabelList& faceCells() const
			{
				return fvPatch::faceCells();
			}

			//- Return delta (P to N) vectors across coupled patch
			virtual tmp<vectorField> delta() const = 0;


		// Interface transfer functions

			//- Return the values of the given internal data adjacent to
			//  the interface as a field
			virtual tmp<labelField> interfaceInternalField
			(
				const unallocLabelList& internalData
			) const = 0;

			//- Initialise interface data transfer
			virtual void initTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& interfaceData
			) const
			{}

			//- Transfer and return neighbour field
			virtual tmp<labelField> transfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& interfaceData
			) const = 0;

			//- Initialise neighbour field transfer
			virtual void initInternalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				unallocLabelList& iF
			) const
			{}

			//- Return neighbour field
			virtual tmp<labelField> internalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& iF
			) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
