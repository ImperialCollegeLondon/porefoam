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
	Foam::cell

Description
	A cell is defined as a list of faces with extra functionality

SourceFiles
	cellI.H

\*---------------------------------------------------------------------------*/

#ifndef cell_H
#define cell_H

#include "faceList.H"
#include "oppositeFace.H"
#include "debugSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class cell;
bool operator==(const cell&, const cell&);
inline bool operator!=(const cell&, const cell&);



class cell
:
	public labelList
{

public:

	// Static data members

		static const char* const typeName;

		//- Debug switch
		static debug::debugSwitch debug;


	// Constructors

		//- Construct null
		inline cell();

		//- Construct given size
		explicit inline cell(label);

		//- Construct from list of labels
		explicit inline cell(const UList<label>&);

		//- Construct by transferring the parameter contents
		explicit inline cell(const Xfer<labelList>&);

		//- Construct from Istream
		inline cell(Istream&);


	// Member Functions

		//- Return number of faces
		inline label nFaces() const;

		//- Return labels of cell vertices
		labelList labels(const unallocFaceList&) const;

		//- Return the cell vertices
		pointField points(const unallocFaceList&, const pointField&) const;

		//- Return cell edges
		edgeList edges(const unallocFaceList&) const;

		//- Return index of opposite face
		label opposingFaceLabel
		(
			const label masterFaceLabel,
			const unallocFaceList& meshFaces
		) const;

		//- Return opposite face oriented the same way as the master face
		oppositeFace opposingFace
		(
			const label masterFaceLabel,
			const unallocFaceList& meshFaces
		) const;


		// Instructions for centre and mag:
		// When one wants to access the cell centre and magnitude, the
		// functionality on the mesh level should be used in preference to the
		// functions provided here. They do not rely to the functionality
		// implemented here, provide additional checking and are more efficient.
		// The cell::centre and cell::mag functions may be removed in the
		// future.

		//- Returns cell centre
		point centre(const pointField&, const unallocFaceList&) const;

		//- Returns cell volume
		scalar mag(const pointField&, const unallocFaceList&) const;


	// Friend Operators

		friend bool operator==(const cell&, const cell&);
		friend bool operator!=(const cell&, const cell&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
