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
	Foam::Function1Types::TableBase

Description
	Base class for table with bounds handling, interpolation and integration

SourceFiles
	TableBase.C

\*---------------------------------------------------------------------------*/

#ifndef TableBase_H
#define TableBase_H

#include "Function1.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationWeights;

namespace Function1Types
{


template<class Type>
class TableBase
:
	public Function1<Type>
{
public:

	// Public data types

		//- Enumeration for handling out-of-bound values
		enum boundsHandling
		{
			ERROR,          //!< Exit with a FatalError
			WARN,           //!< Issue warning and clamp value (default)
			CLAMP,          //!< Clamp value to the start/end value
			REPEAT          //!< Treat as a repeating list
		};


protected:

	// Protected data

		//- Table name
		const word name_;

		//- Enumeration for handling out-of-bound values
		const boundsHandling boundsHandling_;

		//- Interpolation type
		const word interpolationScheme_;

		//- Table data
		List<Tuple2<scalar, Type> > table_;

		//- Extracted values
		mutable autoPtr<scalarField> tableSamplesPtr_;

		//- Interpolator method
		mutable autoPtr<interpolationWeights> interpolatorPtr_;

		//- Cached indices and weights
		mutable labelList currentIndices_;

		mutable scalarField currentWeights_;


	// Protected Member Functions

		//- Return (demand driven) interpolator
		const interpolationWeights& interpolator() const;

		//- Disallow default bitwise assignment
		void operator=(const TableBase<Type>&);


public:

	// Constructors

		//- Construct from dictionary - note table is not populated
		TableBase(const word& name, const dictionary& dict);

		//- Copy constructor. Note: steals interpolator, tableSamples
		TableBase(const TableBase<Type>& tbl);


	//- Destructor
	virtual ~TableBase();


	// Member Functions

		//- Return the out-of-bounds handling as a word
		word boundsHandlingToWord(const boundsHandling& bound) const;

		//- Return the out-of-bounds handling as an enumeration
		boundsHandling wordToBoundsHandling(const word& bound) const;

		//- Set the out-of-bounds handling from enum, return previous setting
		boundsHandling outOfBounds(const boundsHandling& bound);

		//- Check the table for size and consistency
		void check() const;

		//- Check minimum table bounds
		bool checkMinBounds(const scalar x, scalar& xDash) const;

		//- Check maximum table bounds
		bool checkMaxBounds(const scalar x, scalar& xDash) const;

		//- Convert time
		virtual void convertTimeBase(const Time& t);

		//- Return Table value
		virtual Type value(const scalar x) const;

		//- Integrate between two (scalar) values
		virtual Type integrate(const scalar x1, const scalar x2) const;

		//- Return the reference values
		virtual tmp<scalarField> x() const;

		//- Return the dependent values
		virtual tmp<Field<Type> > y() const;

		//- Write all table data in dictionary format
		virtual void writeData(Ostream& os) const;

		//- Write keywords only in dictionary format. Used for non-inline
		//  table types
		virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "TableBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
