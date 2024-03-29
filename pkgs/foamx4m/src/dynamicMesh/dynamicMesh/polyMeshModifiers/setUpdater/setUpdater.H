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
	Foam::setUpdater

Description
	Keeps cell/face/vertex sets up-to-date: both the ones loaded and the ones
	on disk.

SourceFiles
	setUpdater.C
	setUpdaterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef setUpdater_H
#define setUpdater_H

#include "polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class setUpdater
:
	public polyMeshModifier
{
	// Private Member Functions

		//- Updates all sets
		template<class Type>
		void updateSets(const mapPolyMesh& morphMap) const;

		//- Disallow default bitwise copy construct
		setUpdater(const setUpdater&);

		//- Disallow default bitwise assignment
		void operator=(const setUpdater&);


public:

	//- Runtime type information
	TypeName("setUpdater");


	// Constructors

		//- Construct from dictionary
		setUpdater
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyTopoChanger& mme
		);


	// Destructor

		virtual ~setUpdater();


	// Member Functions

		//- Check for topology change
		virtual bool changeTopology() const;

		//- Insert the layer addition/removal instructions
		//  into the topological change
		virtual void setRefinement(polyTopoChange&) const;

		//- Modify motion points to comply with the topological change
		virtual void modifyMotionPoints(pointField& motionPoints) const;

		//- Force recalculation of locally stored data on topological change
		virtual void updateMesh(const mapPolyMesh&);

		//- Write
		virtual void write(Ostream&) const;

		//- Write dictionary
		virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "setUpdaterTemplates.C"
#endif


#endif

// ************************************************************************* //
