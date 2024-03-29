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
	Foam::pointToCell

Description
	A topoSetSource to select cells based on usage of points.

SourceFiles
	pointToCell.C

\*---------------------------------------------------------------------------*/

#ifndef pointToCell_H
#define pointToCell_H

#include "topoSetSource.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class pointToCell
:
	public topoSetSource
{
public:
		//- Enumeration defining the valid options
		enum pointAction
		{
			ANY,    // Cells using any point in set
			ALL     // Cells whose all points are in set
		};

private:

		//- Add usage string
		static addToUsageTable usage_;

		static const NamedEnum<pointAction, 2> pointActionNames_;

		//- Name of set to use
		word setName_;

		//- Option
		pointAction option_;


	// Private Member Functions

		//- Depending on point-to-cell option add to or delete from cellSet.
		void combine(topoSet& set, const bool add) const;


public:

	//- Runtime type information
	TypeName("pointToCell");

	// Constructors

		//- Construct from components
		pointToCell
		(
			const polyMesh& mesh,
			const word& setName,
			const pointAction option
		);

		//- Construct from dictionary
		pointToCell
		(
			const polyMesh& mesh,
			const dictionary& dict
		);

		//- Construct from Istream
		pointToCell
		(
			const polyMesh& mesh,
			Istream&
		);


	// Destructor

		virtual ~pointToCell();


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
