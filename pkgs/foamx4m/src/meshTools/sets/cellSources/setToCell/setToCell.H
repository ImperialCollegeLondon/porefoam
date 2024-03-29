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
	Foam::setToCell

Description
	A topoSetSource to select faces based on cellSet.

SourceFiles
	setToCell.C

\*---------------------------------------------------------------------------*/

#ifndef setToCell_H
#define setToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class setToCell
:
	public topoSetSource
{

	// Private data

		//- Add usage string
		static addToUsageTable usage_;

		//- Name of cellSet
		word setName_;


	// Private Member Functions

		void combine(topoSet& set, const bool add) const;


public:

	//- Runtime type information
	TypeName("setToCell");

	// Constructors

		//- Construct from components
		setToCell
		(
			const polyMesh& mesh,
			const word& setName
		);

		//- Construct from dictionary
		setToCell
		(
			const polyMesh& mesh,
			const dictionary& dict
		);

		//- Construct from Istream
		setToCell
		(
			const polyMesh& mesh,
			Istream&
		);


	// Destructor

		virtual ~setToCell();


	// Member Functions

		virtual void applyToSet
		(
			const topoSetSource::setAction action,
			topoSet&
		) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
