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
	Foam::memInfo

Description
	Memory usage information for the process running this object.

Note
	Uses the information from /proc/\<pid\>/status

SourceFiles
	memInfo.C

\*---------------------------------------------------------------------------*/

#ifndef memInfo_H
#define memInfo_H

#include "OSspecific.H"
#include "POSIX.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class memInfo;

Istream& operator>>(Istream&, memInfo&);
Ostream& operator<<(Ostream&, const memInfo&);



class memInfo
{
	// Private data

		//- Peak memory used by the process (VmPeak in /proc/\<pid\>/status)
		int peak_;

		//- Memory used by the process (VmSize in /proc/\<pid\>/status)
		int size_;

		//- Resident set size of the process (VmRSS in /proc/\<pid\>/status)
		int rss_;


public:

	// Constructors

		//- Construct null
		memInfo();


	//- Destructor
	~memInfo();


	// Member Functions

		//- Parse /proc/\<pid\>/status
		const memInfo& update();

		// Access

			//- Access the stored peak memory (VmPeak in /proc/\<pid\>/status)
			//  The value is stored from the previous update()
			int peak() const
			{
				return peak_;
			}

			//- Access the stored memory size (VmSize in /proc/\<pid\>/status)
			//  The value is stored from the previous update()
			int size() const
			{
				return size_;
			}

			//- Access the stored rss value (VmRSS in /proc/\<pid\>/status)
			//  The value is stored from the previous update()
			int rss() const
			{
				return rss_;
			}

			//- True if the memory information appears valid
			bool valid() const;


	// IOstream Operators

		//- Read peak/size/rss from stream
		friend Istream& operator>>(Istream&, memInfo&);

		//- Write peak/size/rss to stream
		friend Ostream& operator<<(Ostream&, const memInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
