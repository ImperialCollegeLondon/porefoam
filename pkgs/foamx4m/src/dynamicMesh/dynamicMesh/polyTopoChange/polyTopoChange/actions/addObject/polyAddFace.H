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
	Foam::polyAddFace

Description
	A face addition data class. A face can be inflated either from a
	point or from another face and can either be in internal or a
	boundary face.

\*---------------------------------------------------------------------------*/

#ifndef polyAddFace_H
#define polyAddFace_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "label.H"
#include "face.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyAddFace
:
	public topoAction
{
	// Private data

		//- Face identifier
		face face_;

		//- Face owner
		label owner_;

		//- Face neighbour
		label neighbour_;

		//- Master point ID for faces blown up from points
		label masterPointID_;

		//- Master edge ID for faces blown up from edges
		label masterEdgeID_;

		//- Master face ID for faces blown up from faces
		label masterFaceID_;

		//- Does the face flux need to be flipped
		bool flipFaceFlux_;

		//- Boundary patch ID
		label patchID_;

		//- Face zone ID
		label zoneID_;

		//- Face zone flip
		bool zoneFlip_;


	// Private Member Functions

		//- Check
		void check() const
		{
			if (face_.size() < 3)
			{
				FatalErrorIn("polyAddFace::check()")
					<< "Invalid face: less than 3 points.  "
					<< "This is not allowed.\n"
					<< "Face: " << face_
					<< " masterPointID:" << masterPointID_
					<< " masterEdgeID:" << masterEdgeID_
					<< " masterFaceID:" << masterFaceID_
					<< " patchID:" << patchID_
					<< " owner:" << owner_
					<< " neighbour:" << neighbour_
					<< abort(FatalError);
			}

			if (min(face_) < 0)
			{
				FatalErrorIn("polyAddFace::check()")
					<< "Face contains invalid vertex ID: " << face_ << ".  "
					<< "This is not allowed.\n"
					<< "Face: " << face_
					<< " masterPointID:" << masterPointID_
					<< " masterEdgeID:" << masterEdgeID_
					<< " masterFaceID:" << masterFaceID_
					<< " patchID:" << patchID_
					<< " owner:" << owner_
					<< " neighbour:" << neighbour_
					<< abort(FatalError);
			}

			if (min(owner_, neighbour_) >= 0 && owner_ == neighbour_)
			{
				FatalErrorIn("polyAddFace::check()")
					<< "Face owner and neighbour are identical.  "
					<< "This is not allowed.\n"
					<< "Face: " << face_
					<< " masterPointID:" << masterPointID_
					<< " masterEdgeID:" << masterEdgeID_
					<< " masterFaceID:" << masterFaceID_
					<< " patchID:" << patchID_
					<< " owner:" << owner_
					<< " neighbour:" << neighbour_
					<< abort(FatalError);
			}

			if (neighbour_ >= 0 && patchID_ >= 0)
			{
				FatalErrorIn("polyAddFace::check()")
					<< "Patch face has got a neighbour.  Patch ID: " << patchID_
					<< ".  This is not allowed.\n"
					<< "Face: " << face_
					<< " masterPointID:" << masterPointID_
					<< " masterEdgeID:" << masterEdgeID_
					<< " masterFaceID:" << masterFaceID_
					<< " patchID:" << patchID_
					<< " owner:" << owner_
					<< " neighbour:" << neighbour_
					<< abort(FatalError);
			}

			if (owner_ < 0 && zoneID_ < 0)
			{
				FatalErrorIn("polyAddFace::check()")
					<< "Face has no owner and is not in a zone.  "
					<< "This is not allowed.\n"
					<< "Face: " << face_
					<< "Face: " << face_
					<< " masterPointID:" << masterPointID_
					<< " masterEdgeID:" << masterEdgeID_
					<< " masterFaceID:" << masterFaceID_
					<< " patchID:" << patchID_
					<< " owner:" << owner_
					<< " neighbour:" << neighbour_
					<< abort(FatalError);
			}

			if (zoneID_ == -1 && zoneFlip_)
			{
				FatalErrorIn("polyAddFace::check()")
					<< "Specified zone flip for a face that does not  "
					<< "belong to zone.  This is not allowed.\n"
					<< "Face: " << face_
					<< " masterPointID:" << masterPointID_
					<< " masterEdgeID:" << masterEdgeID_
					<< " masterFaceID:" << masterFaceID_
					<< " patchID:" << patchID_
					<< " owner:" << owner_
					<< " neighbour:" << neighbour_
					<< abort(FatalError);
			}
		}


public:

	// Static data members

		//- Runtime type information
		TypeName("addFace");


	// Constructors

		//- Construct null.  Used for constructing lists
		polyAddFace()
		:
			face_(0),
			owner_(-1),
			neighbour_(-1),
			masterPointID_(-1),
			masterEdgeID_(-1),
			masterFaceID_(-1),
			flipFaceFlux_(false),
			patchID_(-1),
			zoneID_(-1),
			zoneFlip_(false)
		{}


		//- Construct from components
		polyAddFace
		(
			const face& f,
			const label owner,
			const label neighbour,
			const label masterPointID,
			const label masterEdgeID,
			const label masterFaceID,
			const bool flipFaceFlux,
			const label patchID,
			const label zoneID,
			const bool zoneFlip
		)
		:
			face_(f),
			owner_(owner),
			neighbour_(neighbour),
			masterPointID_(masterPointID),
			masterEdgeID_(masterEdgeID),
			masterFaceID_(masterFaceID),
			flipFaceFlux_(flipFaceFlux),
			patchID_(patchID),
			zoneID_(zoneID),
			zoneFlip_(zoneFlip)
		{
			check();
		}

		//- Construct from Istream
		polyAddFace(Istream& is)
		:
			face_(is),
			owner_(readLabel(is)),
			neighbour_(readLabel(is)),
			masterPointID_(readLabel(is)),
			masterEdgeID_(readLabel(is)),
			masterFaceID_(readLabel(is)),
			flipFaceFlux_(readBool(is)),
			patchID_(readLabel(is)),
			zoneID_(readLabel(is)),
			zoneFlip_(readBool(is))
		{
			check();
		}

		//- Construct and return a clone
		virtual autoPtr<topoAction> clone() const
		{
			return autoPtr<topoAction>(new polyAddFace(*this));
		}


	// Default Destructor

	// Member Functions

		//- Return face
		const face& newFace() const
		{
			return face_;
		}

		//- Return owner cell
		label owner() const
		{
			return owner_;
		}

		//- Return neighour cell
		label neighbour() const
		{
			return neighbour_;
		}

		//- Is the face mastered by a point
		bool isPointMaster() const
		{
			return masterPointID_ >= 0;
		}

		//- Is the face mastered by an edge
		bool isEdgeMaster() const
		{
			return masterEdgeID_ >= 0;
		}

		//- Is the face mastered by another face
		bool isFaceMaster() const
		{
			return masterFaceID_ >= 0;
		}

		//- Is the face appended with no master
		bool appended() const
		{
			return !isPointMaster() && !isEdgeMaster() && !isFaceMaster();
		}

		//- Return master point ID
		label masterPointID() const
		{
			return masterPointID_;
		}

		//- Return master edge ID
		label masterEdgeID() const
		{
			return masterEdgeID_;
		}

		//- Return master face ID
		label masterFaceID() const
		{
			return masterFaceID_;
		}

		//- Does the face flux need to be flipped
		bool flipFaceFlux() const
		{
			return flipFaceFlux_;
		}

		//- Does the face belong to a boundary patch?
		bool isInPatch() const
		{
			return patchID_ >= 0;
		}

		//- Boundary patch ID
		label patchID() const
		{
			return patchID_;
		}

		//- Does the face belong to a zone?
		bool isInZone() const
		{
			return zoneID_ >= 0;
		}

		//- Is the face only a zone face (i.e. not belonging to a cell)
		bool onlyInZone() const
		{
			return zoneID_ >= 0 && owner_ < 0 && neighbour_ < 0;
		}

		//- Face zone ID
		label zoneID() const
		{
			return zoneID_;
		}

		//- Face zone flip
		label zoneFlip() const
		{
			return zoneFlip_;
		}


	// IOstream Operators

		friend Ostream& operator<<(Ostream& os , const polyAddFace& paf)
		{
			os  << paf.face_
				<< paf.owner_
				<< paf.neighbour_
				<< paf.masterPointID_
				<< paf.masterEdgeID_
				<< paf.masterFaceID_
				<< paf.flipFaceFlux_
				<< paf.patchID_
				<< paf.zoneID_
				<< paf.zoneFlip_;

			return os;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
