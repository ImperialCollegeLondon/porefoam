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

InClass
	Foam::memberFunctionSelectionTables

Description
	Macros for easy insertion into member function selection tables

\*---------------------------------------------------------------------------*/

#ifndef addToMemberFunctionSelectionTable_H
#define addToMemberFunctionSelectionTable_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// add to hash-table of functions with typename as the key
#define addToMemberFunctionSelectionTable\
(baseType,thisType,memberFunction,argNames)                                   \
					                                                          \
	/* Add the thisType memberFunction to the table */                        \
	baseType::add##memberFunction##argNames##MemberFunctionToTable<thisType>  \
	add##thisType##memberFunction##argNames##MemberFunctionTo##baseType##Table_



// add to hash-table of functions with 'lookup' as the key
#define addNamedToMemberFunctionSelectionTable\
(baseType,thisType,memberFunction,argNames,lookup)                            \
					                                                          \
	/* Add the thisType memberFunction to the table, find by lookup name */   \
	baseType::add##memberFunction##argNames##MemberFunctionToTable<thisType>  \
	add_##lookup##_##thisType##memberFunction##argNames##MemberFunctionTo##baseType##Table_(#lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// add to hash-table of functions with typename as the key
// use when baseType doesn't need a template argument (eg, is a typedef)
#define addTemplateToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames)                              \
					                                                          \
	/* Add the thisType memberFunction to the table */                        \
	baseType::add##memberFunction##argNames##MemberFunctionToTable<thisType<Targ> > \
	add##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Table_


// add to hash-table of functions with 'lookup' as the key
// use when baseType doesn't need a template argument (eg, is a typedef)
#define addNamedTemplateToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames,lookup)                       \
					                                                          \
	/* Add the thisType memberFunction to the table, find by lookup name */   \
	baseType::add##memberFunction##argNames##MemberFunctionToTable<thisType<Targ> > \
	add_##lookup##_##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Table_(#lookup)

// use when baseType requires the Targ template argument as well
#define addTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames)                              \
					                                                          \
	/* Add the thisType memberFunction to the table */                        \
	baseType<Targ>::add##memberFunction##argNames##MemberFunctionToTable<thisType<Targ> > \
	add##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Targ##Table_

// use when baseType requires the Targ template argument as well
#define addNamedTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames,lookup)                       \
					                                                          \
	/* Add the thisType memberFunction to the table, find by lookup name */   \
	baseType<Targ>::add##memberFunction##argNames##MemberFunctionToTable<thisType<Targ> > \
	add_##lookup##_##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Targ##Table_(#lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// add to hash-table of functions with typename as the key
// use when baseType requires the Targ template argument as well
#define addTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames)                              \
					                                                          \
	/* Add the thisType memberFunction to the table */                        \
	baseType<Targ>::add##memberFunction##argNames##MemberFunctionToTable<thisType<Targ> > \
	add##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Targ##Table_


// add to hash-table of functions with 'lookup' as the key
// use when baseType requires the Targ template argument as well
#define addNamedTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames,lookup)                       \
					                                                          \
	/* Add the thisType memberFunction to the table, find by lookup name */   \
	baseType<Targ>::add##memberFunction##argNames##MemberFunctionToTable<thisType<Targ> > \
	add_##lookup##_##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Targ##Table_(#lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
