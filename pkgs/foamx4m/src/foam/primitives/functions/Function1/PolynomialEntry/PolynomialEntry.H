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
	Foam::Function1Types::PolynomialEntry

Description
	PolynomialEntry container data entry for scalars. Items are stored in a
	list of Tuple2's. Data is input in the form,
	e.g. for an entry \<entryName\> that describes y = x^2 + 2x^3

	\verbatim
		<entryName>  polynomial
		(
			(1    2)
			(2    3)
		);
	\endverbatim

SourceFiles
	PolynomialEntry.C

\*---------------------------------------------------------------------------*/

#ifndef PolynomialEntry_H
#define PolynomialEntry_H

#include "Function1.H"
#include "Tuple2.H"
#include "Function1Fwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{


template<class Type>
class Polynomial
:
	public Function1<Type>
{
	// Private data

		//- Polynomial coefficients - list of prefactor, exponent
		List<Tuple2<Type, Type> > coeffs_;

		//- Flag to indicate whether poly can be integrated
		bool canIntegrate_;


	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const Polynomial<Type>&);


public:

	//- Runtime type information
	TypeName("polynomial");


	// Constructors

		Polynomial(const word& entryName, const dictionary& dict);

		//- Construct from components
		Polynomial
		(
			const word& entryName,
			const List<Tuple2<Type, Type> >&
		);

		//- Copy constructor
		Polynomial(const Polynomial& poly);

		//- Construct and return a clone
		virtual tmp<Function1<Type> > clone() const
		{
			return tmp<Function1<Type> >(new Polynomial(*this));
		}


	//- Destructor
	virtual ~Polynomial();


	// Member Functions

		// Manipulation

			//- Convert time
			virtual void convertTimeBase(const Time& t);


		// Evaluation

			//- Return Polynomial value
			Type value(const scalar x) const;

			//- Integrate between two (scalar) values
			Type integrate(const scalar x1, const scalar x2) const;


		//- Write in dictionary format
		virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "PolynomialEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
