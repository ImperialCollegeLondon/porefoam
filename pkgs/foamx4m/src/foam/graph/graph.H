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
	Foam::graph

Description
	Class to create, store and output qgraph files.

SourceFiles
	graph.C

\*---------------------------------------------------------------------------*/

#ifndef graph_H
#define graph_H

#include "foamString.H"
#include "point.H"
#include "HashPtrTable.H"
#include "curve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class graph;

Ostream& operator<<(Ostream&, const graph&);



class graph
:
	public HashPtrTable<curve>
{
	// Private Data

		string title_;
		string xName_;
		string yName_;

		scalarField x_;


		struct xy
		{
			scalar x_, y_;

			xy()
			{}

		// Friend Operators

			friend bool operator==(const xy& a, const xy& b)
			{
				return equal(a.x_, b.x_) && equal(a.y_, b.y_);
			}

			friend bool operator!=(const xy& a, const xy& b)
			{
				return !(a == b);
			}

			friend Istream& operator>>(Istream& is, xy& xyd)
			{
				is >> xyd.x_ >> xyd.y_;
				return is;
			}

			friend Ostream& operator<<(Ostream& os, const xy& xyd)
			{
				os << xyd.x_ << ' ' << xyd.y_;
				return os;
			}
		};


	// Private Member Functions

		void readCurves(Istream&);


public:

	// Constructors

		//- Construct from title and labels (no curves)
		graph
		(
			const string& title,
			const string& xName,
			const string& yName,
			const scalarField& x
		);

		//- Construct from title, labels and y data for 1 curve
		graph
		(
			const string& title,
			const string& xName,
			const string& yName,
			const scalarField& x,
			const scalarField& y
		);

		//- Construct from Istream given title and labels
		graph
		(
			const string& title,
			const string& xName,
			const string& yName,
			Istream& is
		);

		//- Construct from Istream
		graph(Istream& is);


	// Member functions

		// Access

			const string& title() const
			{
				return title_;
			}

			const string& xName() const
			{
				return xName_;
			}

			const string& yName() const
			{
				return yName_;
			}


			const scalarField& x() const
			{
				return x_;
			}

			scalarField& x()
			{
				return x_;
			}


			const scalarField& y() const;

			scalarField& y();


		// Write

			//- Abstract base class for a graph writer
			class writer
			{

			protected:

				void writeXY
				(
					const scalarField& x,
					const scalarField& y,
					Ostream&
				) const;

			public:

				//- Runtime type information
				TypeName("writer");

				//- Declare run-time constructor selection table
				declareRunTimeSelectionTable
				(
					autoPtr,
					writer,
					word,
					(),
					()
				);

				// Selectors

					//- Return a reference to the selected writer
					static autoPtr<writer> New
					(
					    const word& writeFormat
					);


				// Constructors

					//- Construct null
					writer()
					{}


				// Destructor

					virtual ~writer()
					{}


				// Member Functions

					// Access

					    //- Return the appropriate fileName extension
					    //  for this graph format
					    virtual const word& ext() const = 0;


					// Write

					    //- Write graph in appropriate format
					    virtual void write(const graph&, Ostream&) const = 0;
			};

			//- Write out graph data as a simple table
			void writeTable(Ostream&) const;

			//- Write graph to stream in given format
			void write(Ostream&, const word& format) const;

			//- Write graph to file in given format
			void write(const fileName& fName, const word& format) const;


	// Friend operators

		//- Ostream Operator
		friend Ostream& operator<<(Ostream&, const graph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
