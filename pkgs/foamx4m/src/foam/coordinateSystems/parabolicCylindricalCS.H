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
	Foam::parabolicCylindricalCS

Description
	Parabolic cylindrical coordinate system.

	Notation: u = a.x() v = a.y() z = a.z();

Note
	The maintenance of this class may lag that of the main types.

SourceFiles
	parabolicCylindricalCS.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicCylindricalCS_H
#define parabolicCylindricalCS_H

#include "coordinateSystem.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class parabolicCylindricalCS
:
	public coordinateSystem
{
	// Private data members

		//- Are angles in degrees? (default = true)
		Switch inDegrees_;


protected:

	// Protected Member Functions

		//- Convert from local coordinate system to the global Cartesian system
		//  with optional translation for the origin
		virtual vector localToGlobal(const vector&, bool translate) const;

		//- Convert from local coordinate system to the global Cartesian system
		//  with optional translation for the origin
		virtual tmp<vectorField> localToGlobal
		(
			const vectorField&,
			bool translate
		) const;

		//- Convert from global Cartesian system to the local coordinate system
		//  with optional translation for the origin
		virtual vector globalToLocal(const vector&, bool translate) const;

		//- Convert from global Cartesian system to the local coordinate system
		//  with optional translation for the origin
		virtual tmp<vectorField> globalToLocal
		(
			const vectorField&,
			bool translate
		) const;

public:

	//- Runtime type information
	TypeName("parabolicCylindrical");


	// Constructors

		//- Construct null
		parabolicCylindricalCS(const bool inDegrees = true);

		//- Construct from origin and rotation
		parabolicCylindricalCS
		(
			const word& name,
			const point& origin,
			const coordinateRotation&,
			const bool inDegrees = true
		);

		//- Construct from dictionary
		parabolicCylindricalCS(const word&, const dictionary&);

		//- Return clone
		virtual autoPtr<coordinateSystem> clone() const
		{
			return autoPtr<coordinateSystem>(new parabolicCylindricalCS(*this));
		}


	// Member Functions

		// Global information about the coordinate system

			//- Directions in which the span is limited
			virtual spanInfo spanLimited() const;

			//- Span bounds
			virtual boundBox spanBounds() const;


		//- Are angles in degrees?
		bool inDegrees() const;

		//- Non-const access to inDegrees
		Switch& inDegrees();


		// Write

			//- Write
			virtual void write(Ostream&) const;

			//- Write dictionary
			virtual void writeDict(Ostream&, bool subDict = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
