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
	Foam::lduInterface

Description
	An abstract base class for implicitly-coupled interfaces
	e.g. processor and cyclic patches.

SourceFiles
	lduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef lduInterface_H
#define lduInterface_H

#include "labelField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


class crMatrix;

template<class T>
class autoPtr;


class lduInterface
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		lduInterface(const lduInterface&);

		//- Disallow default bitwise assignment
		void operator=(const lduInterface&);


public:

	//- Runtime type information
	TypeName("lduInterface");


	// Constructors

		//- Construct null
		lduInterface()
		{}


	//- Destructor
	virtual ~lduInterface();


	// Member Functions

		// Access

			//- Return true if interface is coupled
			virtual bool coupled() const = 0;

			//- Return faceCell addressing
			virtual bool master() const = 0;

			//- Return faceCell addressing
			virtual const unallocLabelList& faceCells() const = 0;


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

			//- Initialise transfer of internal field adjacent to the interface
			virtual void initInternalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& iF
			) const
			{}

			//- Transfer and return internal field adjacent to the interface
			virtual tmp<labelField> internalFieldTransfer
			(
				const Pstream::commsTypes commsType,
				const unallocLabelList& iF
			) const = 0;

			//- Initialise transfer of prolongation matrix to the interface
			virtual void initProlongationTransfer
			(
				const Pstream::commsTypes commsType,
				const crMatrix& filteredP
			) const
			{
				Pout<< "lduInterface::initProlongationTransfer()" << endl;
			}

			//- Transfer and return prolongation matrix adjacent to
			//  the interface
			virtual autoPtr<crMatrix> prolongationTransfer
			(
				const Pstream::commsTypes commsType,
				const crMatrix& filteredP
			) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
