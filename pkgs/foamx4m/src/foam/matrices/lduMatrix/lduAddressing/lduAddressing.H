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
	Foam::lduAddressing

Description
	The class contains the addressing required by the lduMatrix: upper, lower
	and losort.

	The addressing can be created in two ways: either with references to
	upper and lower in which case it stores references or from labelLists,
	in which case it stores the addressing itself. Additionally, the losort
	addressing belongs to the class is as on lazy evaluation.

	The ordering of owner addresses is such that the labels are in
	increasing order, with groups of identical labels for edges "owned" by
	the same point. The neighbour labels are also ordered in ascending
	order but only for groups of edges belonging to each point. An example
	is given below:
	@verbatim
		owner     eighbour
		0         1
		0         20
		1         2
		1         21
		2         3
		2         22
		3         4
		3         23
		4         5
		4         24
		5         6
		5         25
		6         7
		6         26
		7         8
		7         27
		8         9
		8         28
		9         10
		9         29
	@endverbatim

	There exists an alternative way of addressing the owner
	list: instead of repeating the same label in the owner list, it is
	possible to address the start of each point neighbours in the
	neighbour list. This reduces the size of owner addressing from a list
	over all edges to a list over all points + 1:

	@verbatim
		Owner start list: 0 2 4 6 8 10 12 14 16 18
	@endverbatim

	We shall use the second form of the addressing for fast lookup
	of edge label from the known owner and neighbour, using the following
	algorithm:
	-# take the owner label and position the start of lookup
	   using the owner start list
	-# loop through all neighbours of this owner (ending at the start of
	  lookup of owner + 1) until the match with current neighbour is found.
	  The index used on the neighbour list for the match is the edge index.

	While owner start addressing allows us to find the edge owned by the
	points, it is also necessary to find the edges for which the point is
	a neighbour. Losort addressing lists the edges neighboured by the
	point and we shall use the same trick as above to address into this
	list. Thus, for every point the losort start gives the address of the
	first face to neighbour this point.

SourceFiles
	lduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef lduAddressing_H
#define lduAddressing_H

#include "labelList.H"
#include "lduSchedule.H"
#include "lduInterfaceFieldPtrsList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class extendedLduAddressing;


class lduAddressing
{
	// Private data

		//- Number of equations
		label size_;


	//- Demand-driven data

		//- Losort addressing
		mutable labelList* losortPtr_;

		//- Owner start addressing
		mutable labelList* ownerStartPtr_;

		//- Losort start addressing
		mutable labelList* losortStartPtr_;


		// Demand-driven data for ILU precondition with p-order fill in (ILUCp)

			//- Extended addressing with p-order fill-in
			mutable PtrList<extendedLduAddressing> extendedAddr_;

		// Demand-driven data for parallelisation of ILU0 and Cholesky
		// preconditioners where the order of factorisation is important

			//- List of coefficients for internal equations (e.g. faces where
			//  at least one cell is internal cell)
			mutable labelList* internalEqnCoeffsPtr_;

			//- List of coefficients where the owner cell is at the boundary and
			//  the neighbour cell is internal. Needed for correct ordering
			mutable labelList* flippedInternalEqnCoeffsPtr_;

			//- PtrList of coefficients for boundary equations for each coupled
			//  interface (e.g. faces where both cells are boundary cells on a
			//  given processor patch)
			mutable PtrList<dynamicLabelList> boundaryEqnCoeffs_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		lduAddressing(const lduAddressing&);

		//- Disallow default bitwise assignment
		void operator=(const lduAddressing&);

		//- Calculate losort
		void calcLosort() const;

		//- Calculate owner start
		void calcOwnerStart() const;

		//- Calculate losort start
		void calcLosortStart() const;

		//- Calculate internal/boundary equation coefficients given a list of
		//  interfaces
		void calcInternalBoundaryEqnCoeffs
		(
			const lduInterfaceFieldPtrsList& lduInterfaces
		) const;


public:

	// Constructor
	lduAddressing(const label nEqns);


	//- Destructor
	virtual ~lduAddressing();


	// Member Functions

		//- Return number of equations
		label size() const
		{
			return size_;
		}

		//- Return number of interfaces
		virtual label nPatches() const = 0;

		//- Return lower addressing
		virtual const unallocLabelList& lowerAddr() const = 0;

		//- Return upper addressing
		virtual const unallocLabelList& upperAddr() const = 0;

		//- Return patch to internal addressing given patch number
		virtual const unallocLabelList& patchAddr
		(
			const label patchNo
		) const = 0;

		// Return patch field evaluation schedule
		virtual const lduSchedule& patchSchedule() const = 0;

		//- Return losort addressing
		const unallocLabelList& losortAddr() const;

		//- Return owner start addressing
		const unallocLabelList& ownerStartAddr() const;

		//- Return losort start addressing
		const unallocLabelList& losortStartAddr() const;

		//- Return off-diagonal index given owner and neighbour label
		label triIndex(const label a, const label b) const;

		//- Return extended addressing given p index
		const extendedLduAddressing& extendedAddr(const label p) const;

		//- Return list of coefficients for internal equations
		const unallocLabelList& internalEqnCoeffs
		(
			const lduInterfaceFieldPtrsList& lduInterfaces
		) const;

		//- Return internal cell/boundary cell face flips
		const unallocLabelList& flippedInternalEqnCoeffs
		(
			const lduInterfaceFieldPtrsList& lduInterfaces
		) const;

		//- Return list of coefficients for boundary equations given list of
		//  ldu interface fields and coupled interface index
		const dynamicLabelList& boundaryEqnCoeffs
		(
			const lduInterfaceFieldPtrsList& lduInterfaces,
			const label intI
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
