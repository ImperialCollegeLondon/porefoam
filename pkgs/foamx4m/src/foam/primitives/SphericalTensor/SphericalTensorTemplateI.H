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

#include "VectorTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
template<class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor()
{}


// Construct given VectorSpace
template<class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor
(
	const VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>& vs
)
:
	VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>(vs)
{}


// Construct given three Cmpts
template<class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor(const Cmpt& stii)
{
	this->v_[II] = stii;
}


// Construct from Istream
template<class Cmpt>
inline SphericalTensor<Cmpt>::SphericalTensor(Istream& is)
:
	VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt&  SphericalTensor<Cmpt>::ii() const
{
	return this->v_[II];
}


template<class Cmpt>
inline Cmpt& SphericalTensor<Cmpt>::ii()
{
	return this->v_[II];
}


template<class Cmpt>
inline const SphericalTensor<Cmpt>& SphericalTensor<Cmpt>::T() const
{
	return *this;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two spherical tensors
template<class Cmpt>
inline SphericalTensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& st2)
{
	return SphericalTensor<Cmpt>(st1.ii()*st2.ii());
}


//- Inner-product between a spherical tensor and a vector
template<class Cmpt>
inline Vector<Cmpt>
operator&(const SphericalTensor<Cmpt>& st, const Vector<Cmpt>& v)
{
	return Vector<Cmpt>
	(
		st.ii()*v.x(),
					   st.ii()*v.y(),
					                  st.ii()*v.z()
	);
}


//- Inner-product between a vector and a spherical tensor
template<class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const SphericalTensor<Cmpt>& st)
{
	return Vector<Cmpt>
	(
		v.x()*st.ii(),
					   v.y()*st.ii(),
					                  v.z()*st.ii()
	);
}


//- Double-dot-product between a spherical tensor and a spherical tensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& st2)
{
	return 3*st1.ii()*st2.ii();
}


//- Division of a scalar by a sphericalTensor
template<class Cmpt>
inline SphericalTensor<Cmpt>
operator/(const scalar s, const SphericalTensor<Cmpt>& st)
{
	return SphericalTensor<Cmpt>(s/st.ii());
}


template<class Cmpt>
inline Cmpt magSqr(const SphericalTensor<Cmpt>& st)
{
	return 3*magSqr(st.ii());
}


//- Return the trace of a spherical tensor
template<class Cmpt>
inline Cmpt tr(const SphericalTensor<Cmpt>& st)
{
	return 3*st.ii();
}


//- Return the spherical part of a spherical tensor, i.e. itself
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const SphericalTensor<Cmpt>& st)
{
	return st;
}


//- Return the determinant of a spherical tensor
template<class Cmpt>
inline Cmpt det(const SphericalTensor<Cmpt>& st)
{
	return st.ii()*st.ii()*st.ii();
}


//- Return the inverse of a spherical tensor
template<class Cmpt>
inline SphericalTensor<Cmpt> inv(const SphericalTensor<Cmpt>& st)
{
	return SphericalTensor<Cmpt>(1.0/st.ii());
}


template<class Cmpt>
class outerProduct<SphericalTensor<Cmpt>, Cmpt>
{
public:

	typedef SphericalTensor<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SphericalTensor<Cmpt> >
{
public:

	typedef SphericalTensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, SphericalTensor<Cmpt> >
{
public:

	typedef SphericalTensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, Vector<Cmpt> >
{
public:

	typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SphericalTensor<Cmpt> >
{
public:

	typedef Vector<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
