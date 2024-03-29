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
	Nastran surface reader.

	- Uses the Ansa "$ANSA_NAME" or the Hypermesh "$HMNAME COMP" extensions
	  to obtain patch names.
	- Handles Nastran short and long formats, but not free format.
	- Properly handles the Nastran compact floating point notation: \n
	@verbatim
		GRID          28        10.20269-.030265-2.358-8
	@endverbatim


\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Do weird things to extract number
static scalar parseNASCoord(const string& s)
{
	size_t expSign = s.find_last_of("+-");

	if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
	{
		scalar mantissa = readScalar(IStringStream(s.substr(0, expSign))());
		scalar exponent = readScalar(IStringStream(s.substr(expSign+1))());

		if (s[expSign] == '-')
		{
			exponent = -exponent;
		}
		return mantissa*pow(10, exponent);
	}
	else
	{
		return readScalar(IStringStream(s)());
	}
}


bool triSurface::readNAS(const fileName& fName)
{
	IFstream is(fName);

	if (!is.good())
	{
		FatalErrorIn("triSurface::readNAS(const fileName&)")
			<< "Cannot read file " << fName
			<< exit(FatalError);
	}

	// coordinates of point
	DynamicList<point> points;
	// Nastran index of point
	dynamicLabelList indices;
	// Faces in terms of Nastran point indices
	DynamicList<labelledTri> faces;
	// From face group to patch
	Map<label> groupToPatch;
	label nPatches = 0;
	// Name for face group
	Map<word> groupToName;

	// Ansa tags. Denoted by $ANSA_NAME. These will appear just before the
	// first use of a type. We read them and store the pshell types which
	// are used to name the patches.
	label ansaId = -1;
	word ansaType;
	string ansaName;

	// A single warning per unrecognized command
	HashSet<word> unhandledCmd;

	while (is.good())
	{
		string line;
		is.getLine(line);

		// Ansa extension
		if (line.substr(0, 10) == "$ANSA_NAME")
		{
			string::size_type sem0 = line.find (';', 0);
			string::size_type sem1 = line.find (';', sem0+1);
			string::size_type sem2 = line.find (';', sem1+1);

			if
			(
				sem0 != string::npos
			 && sem1 != string::npos
			 && sem2 != string::npos
			)
			{
				ansaId = readLabel
				(
					IStringStream(line.substr(sem0+1, sem1-sem0-1))()
				);
				ansaType = line.substr(sem1+1, sem2-sem1-1);

				string nameString;
				is.getLine(ansaName);
				if (ansaName[ansaName.size()-1] == '\r')
				{
					ansaName = ansaName.substr(1, ansaName.size()-2);
				}
				else
				{
					ansaName = ansaName.substr(1, ansaName.size()-1);
				}

				// Info<< "ANSA tag for NastranID:" << ansaId
				//     << " of type " << ansaType
				//     << " name " << ansaName << endl;
			}
		}


		// Hypermesh extension
		// $HMNAME COMP                   1"partName"
		if
		(
			line.substr(0, 12) == "$HMNAME COMP"
		 && line.find ('"') != string::npos
		)
		{
			label groupId = readLabel
			(
				IStringStream(line.substr(16, 16))()
			);

			IStringStream lineStream(line.substr(32));

			string rawName;
			lineStream >> rawName;

			groupToName.insert(groupId, string::validate<word>(rawName));
			Info<< "group " << groupId << " => " << rawName << endl;
		}


		if (line.empty() || line[0] == '$')
		{
			// Skip empty or comment
			continue;
		}

		// Check if character 72 is continuation
		if (line.size() > 72 && line[72] == '+')
		{
			line = line.substr(0, 72);

			while (true)
			{
				string buf;
				is.getLine(buf);

				if (buf.size() > 72 && buf[72]=='+')
				{
					line += buf.substr(8, 64);
				}
				else
				{
					line += buf.substr(8, buf.size()-8);
					break;
				}
			}
		}

		// Read first word
		IStringStream lineStream(line);
		word cmd;
		lineStream >> cmd;

		if (cmd == "CTRIA3")
		{
			label groupId = readLabel(IStringStream(line.substr(16,8))());
			label a = readLabel(IStringStream(line.substr(24,8))());
			label b = readLabel(IStringStream(line.substr(32,8))());
			label c = readLabel(IStringStream(line.substr(40,8))());


			// Convert group into patch
			Map<label>::const_iterator iter = groupToPatch.find(groupId);

			label patchI;
			if (iter == groupToPatch.end())
			{
				patchI = nPatches++;
				groupToPatch.insert(groupId, patchI);
				Info<< "patch " << patchI << " => group " << groupId << endl;
			}
			else
			{
				patchI = iter();
			}

			faces.append(labelledTri(a, b, c, patchI));
		}
		else if (cmd == "CQUAD4")
		{
			label groupId = readLabel(IStringStream(line.substr(16,8))());
			label a = readLabel(IStringStream(line.substr(24,8))());
			label b = readLabel(IStringStream(line.substr(32,8))());
			label c = readLabel(IStringStream(line.substr(40,8))());
			label d = readLabel(IStringStream(line.substr(48,8))());

			// Convert group into patch
			Map<label>::const_iterator iter = groupToPatch.find(groupId);

			label patchI;
			if (iter == groupToPatch.end())
			{
				patchI = nPatches++;
				groupToPatch.insert(groupId, patchI);
				Info<< "patch " << patchI << " => group " << groupId << endl;
			}
			else
			{
				patchI = iter();
			}

			faces.append(labelledTri(a, b, c, patchI));
			faces.append(labelledTri(c, d, a, patchI));
		}
		else if (cmd == "PSHELL")
		{
			// Read shell type since group gives patchnames
			label groupId = readLabel(IStringStream(line.substr(8,8))());
			if (groupId == ansaId && ansaType == "PSHELL")
			{
				groupToName.insert(groupId, string::validate<word>(ansaName));
				Info<< "group " << groupId << " => " << ansaName << endl;
			}
		}
		else if (cmd == "GRID")
		{
			label index = readLabel(IStringStream(line.substr(8,8))());
			scalar x = parseNASCoord(line.substr(24, 8));
			scalar y = parseNASCoord(line.substr(32, 8));
			scalar z = parseNASCoord(line.substr(40, 8));

			indices.append(index);
			points.append(point(x, y, z));
		}
		else if (cmd == "GRID*")
		{
			// Long format is on two lines with '*' continuation symbol
			// on start of second line.
			// Typical line (spaces compacted)
			// GRID*      126   0 -5.55999875E+02 -5.68730474E+02
			// *         2.14897901E+02

			label index = readLabel(IStringStream(line.substr(8,16))());
			scalar x = parseNASCoord(line.substr(40, 16));
			scalar y = parseNASCoord(line.substr(56, 16));

			is.getLine(line);
			if (line[0] != '*')
			{
				FatalErrorIn("triSurface::readNAS(const fileName&)")
					<< "Expected continuation symbol '*' when reading GRID*"
					<< " (double precision coordinate) output" << nl
					<< "Read:" << line << nl
					<< "File:" << is.name()
					<< " line:" << is.lineNumber()
					<< exit(FatalError);
			}
			scalar z = parseNASCoord(line.substr(8, 16));

			indices.append(index);
			points.append(point(x, y, z));
		}
		else if (unhandledCmd.insert(cmd))
		{
			Info<< "Unhandled Nastran command " << line << nl
				<< "File:" << is.name() << " line:" << is.lineNumber() << endl;
		}
	}

	points.shrink();
	indices.shrink();
	faces.shrink();


	Info<< "Read triangles:" << faces.size() << " points:" << points.size()
		<< endl;

	{
		// Build inverse mapping (index to point)
		Map<label> indexToPoint(2*indices.size());
		forAll(indices, i)
		{
			indexToPoint.insert(indices[i], i);
		}

		// Relabel faces
		forAll(faces, i)
		{
			labelledTri& f = faces[i];

			f[0] = indexToPoint[f[0]];
			f[1] = indexToPoint[f[1]];
			f[2] = indexToPoint[f[2]];
		}
	}


	// Convert groupToPatch to patchList.
	geometricSurfacePatchList patches(nPatches);

	forAllConstIter(Map<word>, groupToName, iter)
	{
		label patchI = groupToPatch[iter.key()];

		patches[patchI] = geometricSurfacePatch
		(
			"empty",
			iter(),
			patchI
		);
	}

	Info<< "patches:" << patches << endl;

	// Transfer DynamicLists to straight ones.
	pointField allPoints(points.xfer());

	// Create triSurface
	*this = triSurface(faces, patches, allPoints, true);

	return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
