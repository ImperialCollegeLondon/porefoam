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

#include "MeshedSurfaceProxy.H"

#include "foamTime.H"
#include "surfMesh.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::MeshedSurfaceProxy<Face>::writeTypes()
{
	return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


template<class Face>
bool Foam::MeshedSurfaceProxy<Face>::canWriteType
(
	const word& ext,
	const bool verbose
)
{
	return checkSupport(writeTypes(), ext, verbose, "writing");
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
	const fileName& name,
	const MeshedSurfaceProxy& surf
)
{
	if (debug)
	{
		Info<< "MeshedSurfaceProxy::write"
			"(const fileName&, const MeshedSurfaceProxy&) : "
			"writing to " << name
			<< endl;
	}

	word ext = name.ext();

	typename writefileExtensionMemberFunctionTable::iterator mfIter =
		writefileExtensionMemberFunctionTablePtr_->find(ext);

	if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
	{
		FatalErrorIn
		(
			"MeshedSurfaceProxy::write(const fileName&)"
		)   << "Unknown file extension " << ext << nl << nl
			<< "Valid types are :" << endl
			<< writeTypes()
			<< exit(FatalError);
	}

	mfIter()(name, surf);
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
	const Time& t,
	const word& surfName
) const
{
	// the surface name to be used
	word name(surfName.size() ? surfName : surfaceRegistry::defaultName);

	if (debug)
	{
		Info<< "MeshedSurfaceProxy::write"
			"(const Time&, const word&) : "
			"writing to " << name
			<< endl;
	}


	// the local location
	const fileName objectDir
	(
		t.timePath()/surfaceRegistry::prefix/name/surfMesh::meshSubDir
	);

	if (!isDir(objectDir))
	{
		mkDir(objectDir);
	}


	// write surfMesh/points
	{
		pointIOField io
		(
			IOobject
			(
				"points",
				t.timeName(),
				surfMesh::meshSubDir,
				t,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		OFstream os
		(
			objectDir/io.name(),
			ios_base::out|ios_base::trunc,
			t.writeFormat(),
			IOstream::currentVersion,
			t.writeCompression()
		);

		io.writeHeader(os);

		os  << this->points();

		io.writeEndDivider(os);
	}


	// write surfMesh/faces
	{
		faceIOList io
		(
			IOobject
			(
				"faces",
				t.timeName(),
				surfMesh::meshSubDir,
				t,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		OFstream os
		(
			objectDir/io.name(),
			ios_base::out|ios_base::trunc,
			t.writeFormat(),
			IOstream::currentVersion,
			t.writeCompression()
		);
		io.writeHeader(os);

		if (this->useFaceMap())
		{
			// this is really a bit annoying (and wasteful) but no other way
			os  << reorder(this->faceMap(), this->faces());
		}
		else
		{
			os  << this->faces();
		}

		io.writeEndDivider(os);
	}


	// write surfMesh/surfZones
	{
		surfZoneIOList io
		(
			IOobject
			(
				"surfZones",
				t.timeName(),
				surfMesh::meshSubDir,
				t,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		// write as ascii
		OFstream os(objectDir/io.name());
		io.writeHeader(os);

		os  << this->surfZones();

		io.writeEndDivider(os);
	}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::MeshedSurfaceProxy
(
	const pointField& pointLst,
	const List<Face>& faceLst,
	const List<surfZone>& zoneLst,
	const labelList& faceMap
)
:
	points_(pointLst),
	faces_(faceLst),
	zones_(zoneLst),
	faceMap_(faceMap)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::~MeshedSurfaceProxy()
{}


// ************************************************************************* //
