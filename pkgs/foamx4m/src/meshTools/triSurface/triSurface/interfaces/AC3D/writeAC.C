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

Description

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeAC(Ostream& os) const
{
	// Write with patches as separate objects under "world" object.
	// Header is taken over from sample file.
	// Defines separate materials for all patches. Recycle colours.

	// Define 8 standard colours as r,g,b components
	static scalar colourMap[] =
	{
		1, 1, 1,
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
		1, 1, 0,
		0, 1, 1,
		1, 0, 1,
		0.5, 0.5, 1
	};

	// Calculate patch face indexing

	labelList faceMap;

	surfacePatchList myPatches(calcPatches(faceMap));


	// Write header. Define materials.

	os  << "AC3Db" << endl;

	forAll(myPatches, patchI)
	{
		const word& pName = myPatches[patchI].name();

		label colourI = patchI % 8;
		label colourCompI = 3 * colourI;

		os  << "MATERIAL \"" << pName << "Mat\" rgb "
			<< colourMap[colourCompI] << ' ' << colourMap[colourCompI+1]
			<< ' ' << colourMap[colourCompI+2]
			<< "  amb 0.2 0.2 0.2  emis 0 0 0  spec 0.5 0.5 0.5  shi 10"
			<< "  trans 0"
			<< endl;
	}

	os  << "OBJECT world" << endl
		<< "kids " << myPatches.size() << endl;


	// Write patch points & faces.

	label faceIndex = 0;

	forAll(myPatches, patchI)
	{
		const surfacePatch& sp = myPatches[patchI];

		os  << "OBJECT poly" << endl
			<< "name \"" << sp.name() << '"' << endl;

		// Create patch with only patch faces included for ease of addressing

		boolList include(size(), false);

		forAll(sp, patchFaceI)
		{
			const label faceI = faceMap[faceIndex++];

			include[faceI] = true;
		}

		labelList pointMap;
		labelList faceMap;

		triSurface patch = subsetMesh(include, pointMap, faceMap);

		// Now we have triSurface for this patch alone. Write it.

		os << "numvert " << patch.nPoints() << endl;

		forAll(patch.localPoints(), ptI)
		{
			const point& pt = patch.localPoints()[ptI];

			os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
		}

		os << "numsurf " << patch.localFaces().size() << endl;

		forAll(patch.localFaces(), faceI)
		{
			const labelledTri& f = patch.localFaces()[faceI];

			os  << "SURF 0x20" << endl          // polygon
				<< "mat " << patchI << endl
				<< "refs " << f.size() << endl;

			os << f[0] << " 0 0" << endl;
			os << f[1] << " 0 0" << endl;
			os << f[2] << " 0 0" << endl;
		}

		os << "kids 0" << endl;
	}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
