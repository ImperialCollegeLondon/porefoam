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

#include "TimeActivatedExplicitSourceTemplate.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Foam::word& Foam::TimeActivatedExplicitSource<Type>::name() const
{
	return name_;
}


template<class Type>
inline const Foam::fvMesh& Foam::TimeActivatedExplicitSource<Type>::mesh() const
{
	return mesh_;
}


template<class Type>
inline bool Foam::TimeActivatedExplicitSource<Type>::active() const
{
	return active_;
}


template<class Type>
inline Foam::scalar Foam::TimeActivatedExplicitSource<Type>::timeStart() const
{
	return timeStart_;
}


template<class Type>
inline Foam::scalar Foam::TimeActivatedExplicitSource<Type>::duration() const
{
	return duration_;
}


template<class Type>
inline Foam::scalar Foam::TimeActivatedExplicitSource<Type>::timeEnd() const
{
	return timeStart_ + duration_;
}


template<class Type>
inline const typename Foam::TimeActivatedExplicitSource<Type>::volumeModeType&
Foam::TimeActivatedExplicitSource<Type>::volumeMode() const
{
	return volumeMode_;
}


template<class Type>
inline const typename Foam::TimeActivatedExplicitSource<Type>::
selectionModeType&
Foam::TimeActivatedExplicitSource<Type>::selectionMode() const
{
	return selectionMode_;
}


template<class Type>
inline const Foam::List<Foam::point>&
Foam::TimeActivatedExplicitSource<Type>::points() const
{
	return points_;
}


template<class Type>
inline const Foam::word&
Foam::TimeActivatedExplicitSource<Type>::cellSetName() const
{
	return cellSetName_;
}


template<class Type>
inline const Foam::scalarList&
Foam::TimeActivatedExplicitSource<Type>::V() const
{
	return V_;
}


template<class Type>
inline const Foam::cellSet&
Foam::TimeActivatedExplicitSource<Type>::cells() const
{
	return cellsPtr_();
}


template<class Type>
inline const Foam::List<typename Foam::TimeActivatedExplicitSource<Type>::
fieldNameValuePair>&
Foam::TimeActivatedExplicitSource<Type>::fieldData() const
{
	return fieldData_;
}


template<class Type>
inline const Foam::labelList&
Foam::TimeActivatedExplicitSource<Type>::fieldIds() const
{
	return fieldIds_;
}


template<class Type>
inline Foam::word& Foam::TimeActivatedExplicitSource<Type>::name()
{
	return name_;
}


template<class Type>
inline bool& Foam::TimeActivatedExplicitSource<Type>::active()
{
	return active_;
}


template<class Type>
inline Foam::scalar& Foam::TimeActivatedExplicitSource<Type>::timeStart()
{
	return timeStart_;
}


template<class Type>
inline Foam::scalar& Foam::TimeActivatedExplicitSource<Type>::duration()
{
	return duration_;
}


template<class Type>
inline typename Foam::TimeActivatedExplicitSource<Type>::volumeModeType&
Foam::TimeActivatedExplicitSource<Type>::volumeMode()
{
	return volumeMode_;
}


template<class Type>
inline typename Foam::TimeActivatedExplicitSource<Type>::selectionModeType&
Foam::TimeActivatedExplicitSource<Type>::selectionMode()
{
	return selectionMode_;
}


template<class Type>
inline Foam::List<Foam::point>&
Foam::TimeActivatedExplicitSource<Type>::points()
{
	return points_;
}


template<class Type>
inline Foam::word& Foam::TimeActivatedExplicitSource<Type>::cellSetName()
{
	return cellSetName_;
}


template<class Type>
inline Foam::scalarList& Foam::TimeActivatedExplicitSource<Type>::V()
{
	return V_;
}


template<class Type>
inline Foam::cellSet& Foam::TimeActivatedExplicitSource<Type>::cells()
{
	return cellsPtr_();
}


template<class Type>
inline Foam::List<typename Foam::TimeActivatedExplicitSource<Type>::fieldNameValuePair>&
Foam::TimeActivatedExplicitSource<Type>::fieldData()
{
	return fieldData_;
}


template<class Type>
inline Foam::labelList& Foam::TimeActivatedExplicitSource<Type>::fieldIds()
{
	return fieldIds_;
}


// ************************************************************************* //
