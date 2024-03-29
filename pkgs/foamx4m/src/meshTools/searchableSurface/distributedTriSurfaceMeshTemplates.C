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

\*---------------------------------------------------------------------------*/

#include "distributedTriSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Type>
//void Foam::distributedTriSurfaceMesh::getField
//(
//    const word& fieldName,
//    const List<pointIndexHit>& info,
//    List<Type>& values
//) const
//{
//    typedef DimensionedField<Type, triSurfaceGeoMesh> DimensionedSurfField;
//
//
//    // Get query data (= local index of triangle)
//    // ~~~~~~~~~~~~~~
//
//    labelList triangleIndex(info.size());
//    autoPtr<mapDistribute> mapPtr
//    (
//        calcLocalQueries
//        (
//            info,
//            triangleIndex
//        )
//    );
//    const mapDistribute& map = mapPtr();
//
//
//    // Do my tests
//    // ~~~~~~~~~~~
//
//    const DimensionedSurfField& fld = lookupObject<DimensionedSurfField>
//    (
//        fieldName
//    );
//    const triSurface& s = static_cast<const triSurface&>(*this);
//
//    values.setSize(triangleIndex.size());
//
//    forAll(triangleIndex, i)
//    {
//        label triI = triangleIndex[i];
//        values[i] = fld[triI];
//    }
//
//
//    // Send back results
//    // ~~~~~~~~~~~~~~~~~
//
//    map.distribute
//    (
//        Pstream::nonBlocking,
//        List<labelPair>(0),
//        info.size(),
//        map.constructMap(),     // what to send
//        map.subMap(),           // what to receive
//        values
//    );
//}


template<class Type>
void Foam::distributedTriSurfaceMesh::distributeFields
(
	const mapDistribute& map
)
{
	typedef DimensionedField<Type, triSurfaceGeoMesh> DimensionedSurfField;

	HashTable<DimensionedSurfField*> fields
	(
		objectRegistry::lookupClass
		<DimensionedSurfField >()
	);

	forAllIter (typename HashTable<DimensionedSurfField*>, fields, fieldIter)
	{
		DimensionedSurfField& field = *fieldIter();

		label oldSize = field.size();

		map.distribute
		(
			Pstream::nonBlocking,
			List<labelPair>(0),
			map.constructSize(),
			map.subMap(),
			map.constructMap(),
			field
		);

		if (debug)
		{
			Info<< "Mapped " << field.typeName << ' ' << field.name()
				<< " from size " << oldSize << " to size " << field.size()
				<< endl;
		}
	}
}


// ************************************************************************* //
