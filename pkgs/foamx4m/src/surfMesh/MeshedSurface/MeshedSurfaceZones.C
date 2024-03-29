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

#include "MeshedSurface.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::checkZones()
{
	// extra safety, ensure we have at some zones
	// and they cover all the faces - fix start silently
	surfZoneList& zones = this->storedZones();
	if (zones.size())
	{
		label count = 0;
		forAll(zones, zoneI)
		{
			zones[zoneI].start() = count;
			count += zones[zoneI].size();
		}

		if (count < this->size())
		{
			WarningIn
			(
				"MeshedSurface::checkZones()\n"
			)
				<< "more faces " << this->size() << " than zones " << count
				<< " ... extending final zone"
				<< endl;

			zones[zones.size()-1].size() += count - this->size();
		}
		else if (count > this->size())
		{
			FatalErrorIn
			(
				"MeshedSurface::checkZones()\n"
			)
				<< "more zones " << count << " than faces " << this->size()
				<< exit(FatalError);
		}
	}
}


template<class Face>
void Foam::MeshedSurface<Face>::sortFacesAndStore
(
	const Xfer< List<Face> >& unsortedFaces,
	const Xfer< labelList >& zoneIds,
	const bool sorted
)
{
	List<Face>  oldFaces(unsortedFaces);
	labelList zones(zoneIds);

	if (sorted)
	{
		// already sorted - simply transfer faces
		this->storedFaces().transfer(oldFaces);
	}
	else
	{
		// unsorted - determine the sorted order:
		// avoid SortableList since we discard the main list anyhow
		labelList faceMap;
		sortedOrder(zones, faceMap);
		zones.clear();

		// sorted faces
		List<Face> newFaces(faceMap.size());
		forAll(faceMap, faceI)
		{
			// use transfer to recover memory where possible
			newFaces[faceI].transfer(oldFaces[faceMap[faceI]]);
		}
		this->storedFaces().transfer(newFaces);
	}
	zones.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
	const UList<surfZone>& srfZones,
	const bool cullEmpty
)
{
	label nZone = 0;

	surfZoneList& zones = this->storedZones();
	zones.setSize(zones.size());
	forAll(zones, zoneI)
	{
		if (srfZones[zoneI].size() || !cullEmpty)
		{
			zones[nZone] = surfZone(srfZones[zoneI], nZone);
			nZone++;
		}
	}
	zones.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
	const UList<label>& sizes,
	const UList<word>& names,
	const bool cullEmpty
)
{
	label start   = 0;
	label nZone = 0;

	surfZoneList& zones = this->storedZones();
	zones.setSize(sizes.size());
	forAll(zones, zoneI)
	{
		if (sizes[zoneI] || !cullEmpty)
		{
			zones[nZone] = surfZone
			(
				names[zoneI],
				sizes[zoneI],
				start,
				nZone
			);
			start += sizes[zoneI];
			nZone++;
		}
	}
	zones.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
	const UList<label>& sizes,
	const bool cullEmpty
)
{
	label start   = 0;
	label nZone = 0;

	surfZoneList& zones = this->storedZones();
	zones.setSize(sizes.size());
	forAll(zones, zoneI)
	{
		if (sizes[zoneI] || !cullEmpty)
		{
			zones[nZone] = surfZone
			(
				word("zone") + ::Foam::name(nZone),
				sizes[zoneI],
				start,
				nZone
			);
			start += sizes[zoneI];
			nZone++;
		}
	}
	zones.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::removeZones()
{
	this->storedZones().clear();
}


// ************************************************************************* //
