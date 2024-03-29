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
	Foam::ggiSAMGInterface

Description
	AMG selected ggi interface.

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	ggiSAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef ggiSAMGInterface_H
#define ggiSAMGInterface_H

#include "SAMGInterface.H"
#include "ggiLduInterface.H"
#include "mapDistribute.H"
#include "crMatrix.H"
#include "lduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class ggiSAMGInterface
:
	public SAMGInterface,
	virtual public ggiLduInterface
{
	// Private data

		//- Reference tor the ggiLduInterface from which this is
		//  agglomerated
		const ggiLduInterface& fineGgiInterface_;

		//- Zone size
		label zoneSize_;

		//- Zone addressing
		labelList zoneAddressing_;

		//- Processor master faces
		//  Per-processor insertion list of local faces into global zone
		//  The list is created on the master side and passed onto the slave
		//  to allow the slave to insert faces in the same order
		labelListList procMasterFaces_;


		// Parallel communication

			//- Communicator to use for parallel communication
			const label comm_;

			//- Message tag used for sending
			const int tag_;

			//- Map-distribute comms tool
			mutable mapDistribute* mapPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		ggiSAMGInterface(const ggiSAMGInterface&);

		//- Disallow default bitwise assignment
		void operator=(const ggiSAMGInterface&);


		//- Init fast reduce
		void initFastReduce() const;

		//- Return mapDistribute
		const mapDistribute& map() const;


	// Private static data

		//- Processor cluster offset index
		static const long procOffset = 12000000;


public:

	//- Runtime type information
	TypeName("ggi");


	// Constructors

		//- Construct from fine-level interface,
		//  local and neighbour row label
		ggiSAMGInterface
		(
			const lduPrimitiveMesh& lduMesh,
			const crMatrix& interfaceProlongation,
			const lduInterfacePtrsList& coarseInterfaces,
			const lduInterface& fineInterface,
			const crMatrix& nbrInterfaceProlongation
		);


	//- Destructor
	virtual ~ggiSAMGInterface();


	// Member Functions

		// Access

			//- Return true if interface is coupled
			virtual bool coupled() const
			{
				return true;
			}


		// Communications support

			//- Return communicator used for parallel communication
			virtual int comm() const
			{
				return comm_;
			}

			//- Return message tag used for sending
			virtual int tag() const
			{
				return tag_;
			}


		// Coefficient selection

			//- Select the given fine-level coefficients and return
			virtual tmp<scalarField> selectCoeffs
			(
				const scalarField& fineCoeffs
			) const;


		// Interface transfer functions

			//- Fast expand
			//  Note: contains global communications
			//  This function will expand data to zone size with fast comms
			//  filling in only the part of the zone that is required
			//  on the shadow side
			template<class Type>
			tmp<Field<Type> > fastExpand(const UList<Type>&) const;

			//- Fast reduce
			//  Note: contains global communications
			//  This function will expand data to zone size with fast comms
			//  and filter the field to cover the need of a shadow zone
			//  This makes sense on coarse levels because all addressing
			//  is one-on-one (one master, one slave).  On the fine level
			//  each face in the zone may contribute to multiple faces meaning
			// that the data needs to be expanded to zone size
			//  HJ, 24/Jun/2011 and 13/Jun/2016
			template<class Type>
			tmp<Field<Type> > fastReduce(const UList<Type>&) const;

			//- Initialise interface data transfer
			virtual void initTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& interfaceData
			) const;

			//- Transfer and return neighbour field
			virtual tmp<labelField> transfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& interfaceData
			) const;

			//- Initialise neighbour field transfer
			virtual void initInternalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& iF
			) const;

			virtual void initInternalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const scalarField& iF
			) const;

			//- Transfer and return internal field adjacent to the interface
			virtual tmp<labelField> internalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& iF
			) const;

			virtual tmp<scalarField> internalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const scalarField& iF
			) const;

			//- Initialise transfer of prolongation matrix to the interface
			virtual void initProlongationTransfer
			(
				const Pstream::commsTypes commsType,
				const crMatrix& filteredP
			) const;

			//- Transfer and return prolongation matrix adjacent to
			//  the interface
			virtual autoPtr<crMatrix> prolongationTransfer
			(
				const Pstream::commsTypes commsType,
				const crMatrix& filteredP
			) const;


		//- GGI interface functions

			//- Is this the master side?
			virtual bool master() const;

			//- Is this the fine level?
			virtual bool fineLevel() const;

			//- Return shadow patch index
			virtual label shadowIndex() const;

			//- Return shadow interface
			virtual const ggiLduInterface& shadowInterface() const;

			//- Return reference tor the ggiLduInterface from which this is
			//  agglomerated
			const ggiLduInterface& fineGgiInterface() const
			{
				return fineGgiInterface_;
			}

			//- Return interface size
			virtual label interfaceSize() const;

			//- Return zone size
			virtual label zoneSize() const;

			//- Return zone addressing
			virtual const labelList& zoneAddressing() const;

			//- Return addressing
			virtual const labelListList& ggiAddressing() const;

			//- Is the patch localised on a single processor
			virtual bool localParallel() const;

			//- Processor master face insertion list
			const labelListList& procMasterFaces() const;

			//- Return weights
			virtual const scalarListList& ggiWeights() const;

			//- Return face transformation tensor
			virtual const tensorField& forwardT() const;

			//- Return neighbour-cell transformation tensor
			virtual const tensorField& reverseT() const;

			//- Expand addressing to zone
			//  Used in optimised AMG coarsening
			virtual void expandAddrToZone(labelField&) const;

			//- Expand crMatrix to zone
			//  Used in optimised SAMG coarsening
			virtual void expandCrMatrixToZone(crMatrix&) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "ggiSAMGInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
