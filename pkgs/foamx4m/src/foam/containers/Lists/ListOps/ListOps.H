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

InNamspace
	Foam

Description
	Various functions to operate on Lists.

SourceFiles
	ListOps.C
	ListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ListOps_H
#define ListOps_H

#include "labelList.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Renumber the values (not the indices) of a list.
//  Negative ListType elements are left as is.
template<class ListType>
ListType renumber(const UList<label>& oldToNew, const ListType&);

//- Inplace renumber the values of a list.
//  Negative ListType elements are left as is.
template<class ListType>
void inplaceRenumber(const UList<label>& oldToNew, ListType&);


//- Reorder the elements (indices, not values) of a list.
//  Negative ListType elements are left as is.
template<class ListType>
ListType reorder(const UList<label>& oldToNew, const ListType&);

//- Inplace reorder the elements of a list.
//  Negative ListType elements are left as is.
template<class ListType>
void inplaceReorder(const UList<label>& oldToNew, ListType&);


// Variants to work with iterators and sparse tables.
// Need to have iterators and insert()

//- Map values. Do not map negative values.
template<class Container>
void inplaceMapValue(const UList<label>& oldToNew, Container&);

//- Recreate with mapped keys. Do not map elements with negative key.
template<class Container>
void inplaceMapKey(const UList<label>& oldToNew, Container&);


//- Generate the (stable) sort order for the list
template<class T>
void sortedOrder(const UList<T>&, labelList& order);

template<class T, class Cmp>
void sortedOrder(const UList<T>&, labelList& order, const Cmp& cmp);

//- Generate (sorted) indices corresponding to duplicate list values
template<class T>
void duplicateOrder(const UList<T>&, labelList& order);

template<class T, class Cmp>
void duplicateOrder(const UList<T>&, labelList& order, const Cmp& cmp);

//- Generate (sorted) indices corresponding to unique list values
template<class T>
void uniqueOrder(const UList<T>&, labelList& order);

template<class T, class Cmp>
void uniqueOrder(const UList<T>&, labelList& order, const Cmp& cmp);

//- Extract elements of List when select is a certain value.
//  eg, to extract all selected elements:
//    subset<bool, labelList>(selectedElems, true, lst);
template<class T, class ListType>
ListType subset(const UList<T>& select, const T& value, const ListType&);

//- Inplace extract elements of List when select is a certain value.
//  eg, to extract all selected elements:
//    inplaceSubset<bool, labelList>(selectedElems, true, lst);
template<class T, class ListType>
void inplaceSubset(const UList<T>& select, const T& value, ListType&);

//- Extract elements of List when select is true
//  eg, to extract all selected elements:
//    subset<boolList, labelList>(selectedElems, lst);
//  Note a labelHashSet could also be used for the bool-list
template<class BoolListType, class ListType>
ListType subset(const BoolListType& select, const ListType&);

//- Inplace extract elements of List when select is true
//  eg, to extract all selected elements:
//    inplaceSubset<boolList, labelList>(selectedElems, lst);
//  Note a labelHashSet could also be used for the bool-list
template<class BoolListType, class ListType>
void inplaceSubset(const BoolListType& select, ListType&);

//- Invert one-to-one map. Unmapped elements will be -1.
labelList invert(const label len, const UList<label>&);

//- Invert one-to-many map. Unmapped elements will be size 0.
labelListList invertOneToMany(const label len, const UList<label>&);

//- Invert many-to-many.
//  Input and output types need to be inherited from List.
//  eg, faces to pointFaces.
template<class InList, class OutList>
void invertManyToMany(const label len, const UList<InList>&, List<OutList>&);

template<class InList, class OutList>
List<OutList> invertManyToMany(const label len, const UList<InList>& in)
{
	List<OutList> out;
	invertManyToMany<InList,OutList>(len, in, out);
	return out;
}

//- Create identity map (map[i] == i) of given length
labelList identity(const label len);

//- Find first occurence of given element and return index,
//  return -1 if not found. Linear search.
template<class ListType>
label findIndex
(
	const ListType&,
	typename ListType::const_reference,
	const label start=0
);

//- Find all occurences of given element. Linear search.
template<class ListType>
labelList findIndices
(
	const ListType&,
	typename ListType::const_reference,
	const label start=0
);

//- Opposite of findIndices: set values at indices to given value
template<class ListType>
void setValues
(
	ListType&,
	const UList<label>& indices,
	typename ListType::const_reference
);

//- Opposite of findIndices: set values at indices to given value
template<class ListType>
ListType createWithValues
(
	const label sz,
	typename ListType::const_reference initValue,
	const UList<label>& indices,
	typename ListType::const_reference setValue
);

//- Find index of max element (and larger than given element).
//  return -1 if not found. Linear search.
template<class ListType>
label findMax(const ListType&, const label start=0);


//- Find index of min element (and less than given element).
//  return -1 if not found. Linear search.
template<class ListType>
label findMin(const ListType&, const label start=0);


//- Find first occurence of given element in sorted list and return index,
//  return -1 if not found. Binary search.
template<class ListType>
label findSortedIndex
(
	const ListType&,
	typename ListType::const_reference,
	const label start=0
);


//- Find last element < given value in sorted list and return index,
//  return -1 if not found. Binary search.
template<class ListType, class BinaryOp>
label findLower
(
	const ListType&,
	typename ListType::const_reference,
	const label stary,
	const BinaryOp& bop
);


//- Find last element < given value in sorted list and return index,
//  return -1 if not found. Binary search.
template<class ListType>
label findLower
(
	const ListType&,
	typename ListType::const_reference,
	const label start=0
);


//- To construct a List from a C array. Has extra Container type
//  to initialise e.g. wordList from arrays of char*.
template<class Container, class T, int nRows>
List<Container> initList(const T[nRows]);


//- To construct a (square) ListList from a C array. Has extra Container type
//  to initialise e.g. faceList from arrays of labels.
template<class Container, class T, int nRows, int nColumns>
List<Container> initListList(const T[nRows][nColumns]);


//- Helper class for list to append y onto the end of x
template<class T>
class ListAppendEqOp
{
public:
	void operator()(List<T>& x, const List<T>& y) const;
};


//- Reverse a list. First element becomes last element etc.
template<class ListType>
ListType reverseList(const ListType& list);


//- Inplace reversal of a list using Swap.
template<class ListType>
void inplaceReverseList(ListType& list);


//- Rotate a list by n places. If n is positive rotate clockwise/right/down.
//  If n is negative rotate anti-clockwise/left/up.
template<class ListType>
ListType rotateList(const ListType& list, const label n);


//- Inplace reversal of a list using the Reversal Block Swapping algorithm.
template<template<typename> class ListType, class DataType>
void inplaceRotateList(ListType<DataType>& list, label n);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "ListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

