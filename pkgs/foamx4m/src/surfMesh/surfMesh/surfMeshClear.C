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

#include "surfMesh.H"
#include "globalMeshData.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMesh::removeZones()
{
	if (debug)
	{
		Info<< "void surfMesh::removeZones(): "
			<< "Removing surface zones."
			<< endl;
	}

	// Remove the surface zones
	storedZones().clear();

	clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::clearGeom()
{
	if (debug)
	{
		Info<< "void surfMesh::clearGeom() : "
			<< "clearing geometric data"
			<< endl;
	}

	MeshReference::clearGeom();
}


void Foam::surfMesh::clearAddressing()
{
	if (debug)
	{
		Info<< "void surfMesh::clearAddressing() : "
			<< "clearing topology"
			<< endl;
	}

	MeshReference::clearPatchMeshAddr();
}


void Foam::surfMesh::clearOut()
{
	clearGeom();
	clearAddressing();
}


// ************************************************************************* //
