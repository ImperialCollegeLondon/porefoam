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
	Foam::polyMeshModifier

Description
	Virtual base class for mesh modifiers.

SourceFiles
	polyMeshModifier.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshModifier_H
#define polyMeshModifier_H

#include "label.H"
#include "word.H"
#include "Switch.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChanger;
class polyTopoChange;
class mapPolyMesh;

// Forward declaration of friend functions and operators
class polyMeshModifier;

Ostream& operator<<(Ostream&, const polyMeshModifier&);



class polyMeshModifier
{
	// Private data

		//- Name of zone
		word name_;

		//- Index of zone
		label index_;

		//- Reference to morph engine
		const polyTopoChanger& topoChanger_;

		//- Activation switch
		mutable Switch active_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		polyMeshModifier(const polyMeshModifier&);

		//- Disallow default bitwise assignment
		void operator=(const polyMeshModifier&);


public:

	// Static data members

		//- Runtime type information
		TypeName("meshModifier");


	// Declare run-time constructor selection tables

		declareRunTimeSelectionTable
		(
			autoPtr,
			polyMeshModifier,
			dictionary,
			(
				const word& name,
				const dictionary& dict,
				const label index,
				const polyTopoChanger& mme
			),
			(name, dict, index, mme)
		);


	// Constructors

		//- Construct from components
		polyMeshModifier
		(
			const word& name,
			const label index,
			const polyTopoChanger& mme,
			const bool act
		);


	// Selectors

		//- Select constructed from dictionary
		static autoPtr<polyMeshModifier> New
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyTopoChanger& mme
		);


	// Destructor
		virtual ~polyMeshModifier();


	// Member Functions

		//- Return name
		const word& name() const
		{
			return name_;
		}

		//- Return the index of this patch in the boundaryMesh
		label index() const
		{
			return index_;
		}

		//- Return reference to morph engine
		const polyTopoChanger& topoChanger() const;

		//- Check for topology change
		virtual bool changeTopology() const = 0;

		//- Insert the topological change instructions
		virtual void setRefinement(polyTopoChange&) const = 0;

		//- Modify motion points to comply with the topological change
		virtual void modifyMotionPoints(pointField& motionPoints) const = 0;

		//- Force recalculation of locally stored data on topological change
		virtual void updateMesh(const mapPolyMesh&) = 0;


		// Activation and deactivation

			const Switch& active() const
			{
				return active_;
			}

			//- Activate mesh modifier
			void enable() const
			{
				active_ = true;
			}

			//- Activate mesh modifier
			void disable() const
			{
				active_ = false;
			}


		//- Write
		virtual void write(Ostream&) const = 0;

		//- Write dictionary
		virtual void writeDict(Ostream&) const = 0;


	// Ostream Operator

		friend Ostream& operator<<(Ostream&, const polyMeshModifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
