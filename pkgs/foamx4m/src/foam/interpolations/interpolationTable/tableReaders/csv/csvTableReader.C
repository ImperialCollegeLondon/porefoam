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

#include "csvTableReader.H"
#include "IFstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::csvTableReader<Type>::csvTableReader(const dictionary& dict)
:
	tableReader<Type>(dict),
	headerLine_(readBool(dict.lookup("hasHeaderLine"))),
	timeColumn_(readLabel(dict.lookup("timeColumn"))),
	componentColumns_(dict.lookup("valueColumns")),
	separator_(dict.lookupOrDefault<string>("separator", string(","))[0])
{
	if (componentColumns_.size() != pTraits<Type>::nComponents)
	{
		FatalErrorInFunction
			<< componentColumns_ << " does not have the expected length "
			<< pTraits<Type>::nComponents << endl
			<< exit(FatalError);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::csvTableReader<Type>::~csvTableReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
	// doesn't recognize specialization otherwise
	template<>
	scalar csvTableReader<scalar>::readValue(const List<string>& splitted)
	{
		if (componentColumns_[0] >= splitted.size())
		{
			FatalErrorInFunction
				<< "No column " << componentColumns_[0] << " in "
				<< splitted << endl
				<< exit(FatalError);
		}

		return readScalar(IStringStream(splitted[componentColumns_[0]])());
	}


	template<class Type>
	Type csvTableReader<Type>::readValue(const List<string>& splitted)
	{
		Type result;

		for(label i = 0;i < pTraits<Type>::nComponents; i++)
		{
			if (componentColumns_[i] >= splitted.size())
			{
				FatalErrorInFunction
					<< "No column " << componentColumns_[i] << " in "
					<< splitted << endl
					<< exit(FatalError);
			}

			result[i] = readScalar
			(
				IStringStream(splitted[componentColumns_[i]])()
			);
		}

		return result;
	}
}


template<class Type>
void Foam::csvTableReader<Type>::operator()
(
	const fileName& fName,
	List<Tuple2<scalar, Type> >& data
)
{
	IFstream in(fName);

	DynamicList<Tuple2<scalar, Type> > values;

	// Skip header
	if (headerLine_)
	{
		string line;
		in.getLine(line);
	}

	while (in.good())
	{
		string line;
		in.getLine(line);

		DynamicList<string> splitted;

		std::size_t pos = 0;
		while (pos != std::string::npos)
		{
			std::size_t nPos = line.find(separator_, pos);

			if (nPos == std::string::npos)
			{
				splitted.append(line.substr(pos));
				pos=nPos;
			}
			else
			{
				splitted.append(line.substr(pos, nPos-pos));
				pos=nPos+1;
			}
		}

		if (splitted.size() <= 1)
		{
			break;
		}

		scalar time = readScalar(IStringStream(splitted[timeColumn_])());
		Type value = readValue(splitted);

		values.append(Tuple2<scalar,Type>(time, value));
	}

	data.transfer(values);
}


template<class Type>
void Foam::csvTableReader<Type>::operator()
(
	const fileName& fName,
	List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& data
)
{
	NotImplemented;
}


template<class Type>
void Foam::csvTableReader<Type>::write(Ostream& os) const
{
	tableReader<Type>::write(os);

	os.writeKeyword("hasHeaderLine")
		<< headerLine_ << token::END_STATEMENT << nl;
	os.writeKeyword("timeColumn")
		<< timeColumn_ << token::END_STATEMENT << nl;

	// Force writing labelList in ascii
	os.writeKeyword("valueColumns");
	if (os.format() == IOstream::BINARY)
	{
		os.format(IOstream::ASCII);
		os  << componentColumns_;
		os.format(IOstream::BINARY);
	}
	os  << token::END_STATEMENT << nl;

	os.writeKeyword("separator")
		<< string(separator_) << token::END_STATEMENT << nl;
}


// ************************************************************************* //
