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

#include "FTRsurfaceFormat.H"
#include "Keyed.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::FTRsurfaceFormat<Face>::FTRsurfaceFormat
(
	const fileName& filename
)
{
	read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::FTRsurfaceFormat<Face>::read
(
	const fileName& filename
)
{
	this->clear();

	IFstream is(filename);
	if (!is.good())
	{
		FatalErrorIn
		(
			"fileFormats::FTRsurfaceFormat::read(const fileName&)"
		)
			<< "Cannot read file " << filename
			<< exit(FatalError);
	}

	List<ftrPatch> ftrPatches(is);

	// points read directly
	is >> this->storedPoints();

	// triFaces read with attached keys
	List< Keyed<triFace> > facesRead(is);

	List<Face>  faceLst(facesRead.size());
	labelList zoneIds(facesRead.size());

	// disentangle faces/keys - already triangulated
	forAll(facesRead, faceI)
	{
		// unfortunately cannot transfer to save memory
		faceLst[faceI] = facesRead[faceI];
		zoneIds[faceI] = facesRead[faceI].key();
	}

	this->storedFaces().transfer(faceLst);
	this->storedZoneIds().transfer(zoneIds);
	facesRead.clear();

	// change ftrPatch into surfZoneIdentifier
	List<surfZoneIdentifier> newZones(ftrPatches.size());
	forAll(newZones, zoneI)
	{
		newZones[zoneI] = surfZoneIdentifier
		(
			ftrPatches[zoneI].name(),
			zoneI
		);
	}

	this->storedZoneToc().transfer(newZones);
	return true;
}


// ************************************************************************* //
