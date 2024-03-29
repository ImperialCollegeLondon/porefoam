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
	meshObjectBase

Description
	Abstract base-class for dynamic mesh objects used to automate
	their allocation to the mesh database and mesh-related updates

SourceFiles
	meshObjectBase.C

\*---------------------------------------------------------------------------*/

#ifndef meshObjectBase_H
#define meshObjectBase_H

#include "typeInfo.H"
#include "HashTable.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class mapPolyMesh;


class meshObjectBase
{
public:

	//- Runtime type information
	TypeName("meshObject");


	// Static functions

		//- Update topology on all mesh objects
		template<class Mesh>
		static void allUpdateTopology(const Mesh& mesh, const mapPolyMesh& mpm)
		{
			HashTable<const meshObjectBase*> tbl =
				mesh.objectRegistry::template lookupClass<meshObjectBase>();

			if (Mesh::debug)
			{
				InfoIn
				(
					"static void meshObjectBase::"
					"allUpdateTopology(const Mesh& mesh, "
					"const mapPolyMesh& mpm)"
				)   << "Mesh objects to update: " << tbl.toc() << endl;
			}

			for
			(
				HashTable<const meshObjectBase*>::iterator iter =
					tbl.begin();
				iter != tbl.end();
				++iter
			)
			{
				const meshObjectBase& obj = *(iter());

				if (Mesh::debug)
				{
					Info << "Updating object " << obj.type() << endl;
				}

				obj.updateMesh(mpm);
			}

			if (Mesh::debug)
			{
				Info << "Done update topology" << endl;
			}
		}


		//- Move points on all mesh objects
		template<class Mesh>
		static void allMovePoints(const Mesh& mesh)
		{
			if (mesh.moving())
			{
				HashTable<const meshObjectBase*> tbl =
					mesh.objectRegistry::template lookupClass<meshObjectBase>();

				if (Mesh::debug)
				{
					InfoIn
					(
					    "static void meshObjectBase::"
					    "allMovePoints(const Mesh& mesh)"
					)   << "Mesh objects to move: " << tbl.toc() << endl;
				}

				for
				(
					HashTable<const meshObjectBase*>::iterator iter =
					    tbl.begin();
					iter != tbl.end();
					++iter
				)
				{
					const meshObjectBase& obj = *(iter());

					if (Mesh::debug)
					{
					    Info<< "Moving object "
					        << " of type " << obj.type() << endl;
					}

					obj.movePoints();
				}

				if (Mesh::debug)
				{
					Info << "Done moving" << endl;
				}
			}
		}

		//- Move points on all mesh objects
		template<class Mesh>
		static void allDelete(const Mesh& mesh)
		{
			HashTable<const meshObjectBase*> tbl =
				mesh.objectRegistry::template lookupClass<meshObjectBase>();

			for
			(
				HashTable<const meshObjectBase*>::iterator iter =
					tbl.begin();
				iter != tbl.end();
				++iter
			)
			{
				const meshObjectBase& obj = *(iter());

				if (Mesh::debug)
				{
					Info << "Deleting object " << obj.type() << endl;
				}

				obj.deleteObject();
			}
		}


	//- Destructor
	virtual ~meshObjectBase()
	{}


	// Member Functions

		//- Update after mesh motion
		virtual bool movePoints() const = 0;

		//- Update after topology change
		virtual bool updateMesh(const mapPolyMesh&) const = 0;

		//- Delete object from database
		virtual bool deleteObject() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
