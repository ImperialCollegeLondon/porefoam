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
	Foam::curve

Description
	A single curve in a graph.

SourceFiles
	curve.C

\*---------------------------------------------------------------------------*/

#ifndef curve_H
#define curve_H

#include "foamString.H"
#include "primitiveFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class curve;
Ostream& operator<<(Ostream&, const curve&);



class curve
:
	public scalarField
{

public:

	//- The style (line, symbol, etc) of a curve
	class curveStyle
	{

	public:

		//- Enumeration definitions
		enum curveStyleNo
		{
			CONTINUOUS,
			SYMBOL,
			SYMBOL_CURVE,
			SYMBOL_WITH_ERROR_BARS,
			SYMBOL_WITH_VARIABLE_SIZE
		};


	private:

		//- Private data
		curveStyleNo CurveStyleNo;


	public:


		// Constructors

			//- Construct given a curveStyleNo
			curveStyle(const curveStyleNo csn)
			:
				CurveStyleNo(csn)
			{}

			//- Construct from Istream
			curveStyle(Istream& is)
			:
				CurveStyleNo(curveStyleNo(readInt(is)))
			{}


		// Ostream operator

			friend Ostream& operator<<(Ostream& os, const curveStyle& cs)
			{
				os << int(cs.CurveStyleNo);
				return os;
			}
	};


private:

	// private data

		string name_;
		curveStyle style_;


public:

	// Constructors

		//- Construct as interpolation of an existing curve
		//curve(const curve&, const label);

		//- Construct from name, style and size
		curve
		(
			const string& name,
			const curveStyle& style,
			const label l
		);

		//- Construct from the components
		curve
		(
			const string&,
			const curveStyle&,
			const scalarField& y
		);

		//- Create and return a clone
		autoPtr<curve> clone() const
		{
			return autoPtr<curve>(new curve(*this));
		}


	// Member functions

		// Access

			const string& name() const
			{
				return name_;
			}

			const curveStyle& style() const
			{
				return style_;
			}


	// Friend functions

		//- Gradient of the curve
		//friend curve grad(const curve&);


	// Ostream operator

		friend Ostream& operator<<(Ostream&, const curve&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
