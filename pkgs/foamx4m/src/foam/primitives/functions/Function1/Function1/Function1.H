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
	Foam::Function1

Description
	Top level data entry class for use in dictionaries. Provides a mechanism
	to specify a variable as a certain type, e.g. constant or table, and
	provide functions to return the (interpolated) value, and integral between
	limits.

SourceFiles
	Function1.C
	Function1New.C

\*---------------------------------------------------------------------------*/

#ifndef Function1_H
#define Function1_H

#include "dictionary.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

// Forward declaration of friend functions and operators
template<class Type> class Function1;
template<class Type> Ostream& operator<<(Ostream&, const Function1<Type>&);


template<class Type>
class Function1
:
	public refCount
{
	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const Function1<Type>&);


protected:

	// Protected data

		//- Name of entry
		const word name_;


public:

	//- Runtime type information
	TypeName("Function1")

	//- Declare runtime constructor selection table
	declareRunTimeSelectionTable
	(
		autoPtr,
		Function1,
		dictionary,
		(
			const word& entryName,
			const dictionary& dict
		),
		(entryName, dict)
	);


	// Constructor

		//- Construct from entry name
		Function1(const word& entryName);

		//- Copy constructor
		Function1(const Function1<Type>& de);

		//- Construct and return a clone
		virtual tmp<Function1<Type> > clone() const
		{
			return tmp<Function1<Type> >(new Function1<Type>(*this));
		}


	//- Selector
	static autoPtr<Function1<Type> > New
	(
		const word& entryName,
		const dictionary& dict
	);


	//- Destructor
	virtual ~Function1();


	// Member Functions

		// Access

			//- Return the name of the entry
			const word& name() const;


		// Manipulation

			//- Convert time
			virtual void convertTimeBase(const Time& t);


		// Evaluation

			//- Return value as a function of (scalar) independent variable
			virtual Type value(const scalar x) const;

			//- Return value as a function of (scalar) independent variable
			virtual tmp<Field<Type> > value(const scalarField& x) const;

			//- Integrate between two (scalar) values
			virtual Type integrate(const scalar x1, const scalar x2) const;

			//- Integrate between two (scalar) values
			virtual tmp<Field<Type> > integrate
			(
				const scalarField& x1,
				const scalarField& x2
			) const;


		// I/O

			//- Ostream Operator
			friend Ostream& operator<< <Type>
			(
				Ostream& os,
				const Function1<Type>& de
			);

			//- Write in dictionary format
			virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1(Type)                                                    \
					                                                           \
	defineNamedTemplateTypeNameAndDebug(Function1<Type>, 0);                   \
					                                                           \
	defineTemplateRunTimeSelectionTable                                        \
	(                                                                          \
		Function1<Type>,                                                       \
		dictionary                                                             \
	);


#define makeFunction1Type(SS, Type)                                            \
					                                                           \
	defineNamedTemplateTypeNameAndDebug(Function1Types::SS<Type>, 0);          \
					                                                           \
	Function1<Type>::adddictionaryConstructorToTable<Function1Types::SS<Type> > \
		add##SS##Type##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "Function1.C"
	#include "Constant.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
