/*-------------------------------------------------------------------------*\
 Light weight surface utility functions, alternatives to openfoam
 This is part of surfLib, a library for working with surface files and data

 Copyright (C) 2018-2020  Ali Qaseminejad Raeini

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/


#include <vector>
#include "typses.h"
#include "InputFile.h"

using label=int;
using scalar=double;
//using face = std::array<label,4>;
using word = std::string;

class face: public std::array<label,4>
{
 public:
    int zone;
    face() : zone(0) { }
    face(const std::array<label,4>& pointIs) : std::array<label,4> (pointIs) { }
};

template <typename T>  using DynamicList=std::vector<T>;
template <typename T>  using DynamicField=std::vector<T>;
typedef std::vector<std::vector<label> >  labelDynamicListList   ;
typedef std::vector<ints>  labelListList   ;
typedef ints  labelList   ;
typedef piece<face> facePiece    ;
typedef std::array<facePiece,256> facePieceList    ;
typedef dbl3 point    ;
typedef dbl3s vectorField    ;
//typedef vars<point> pointField    ;
struct surfMsh {
	std::vector<face> faces;
	facePieceList faces_bs;
	DynamicField<point> points;
};
#define append  push_back
#define setSize resize
#define Info    std::cout



int  smoothSurf(InputFile& inp, facePieceList& facezs, piece<point> points);
void writeSurfaceFiles(const facePieceList& facezsz, const piece<point> pointsAll, const std::string& fnams);
void writeMergeSurfaceFile(const facePieceList& facezsz, const piece<point> pointsAll, const std::string& fnams);

dbl3s mapFacesGetPoints(std::vector<std::vector<face> >& facess, const piece<point> pointsAll);
std::vector<std::vector<face> > getbsoleteOrderedFaces(const facePieceList& facezsz, int vv);

inline int appendUnique(DynamicList<label>& dynList, label value)
{
	for_i(dynList) if(dynList[i]==value) return 0;
	dynList.append(value);
	return 1;
}

template<typename T, template<typename ...> class C >
 int findPosi(const C<T>& list, const T& val)  {  return distance(list.begin(), find(list.begin(), list.end(),val));  }
template<typename T, template<typename ...> class C>
 int hasValue(const C<T>& list, const T& val)  {  return find(list.begin(), list.end(),val) != list.end(); }
template<class C>
 typename C::const_iterator nextCircIter(const C& f, typename C::const_iterator itr)  {  return ++itr==f.end() ? f.begin() : itr; }
template<class C>
 typename C::const_iterator prevCircIter(const C& f, typename C::const_iterator itr)  {  return itr==f.begin() ? f.end()-1 : --itr; }


inline labelListList getPointPoints(size_t nPoints, const facePiece& faces)
{
    labelDynamicListList pntPntsTmp(nPoints);
    labelListList pntPoints(nPoints);
    for_(faces,fI)
    {
        const face & f=faces[fI];
        for(auto pItr = f.cbegin(); pItr != f.cend(); ++pItr)
        {
           appendUnique(pntPntsTmp[ *pItr ], *nextCircIter(f,pItr));
           appendUnique(pntPntsTmp[ *pItr ], *prevCircIter(f,pItr));
        }
    }

    for_i(pntPoints)
    {
        pntPoints[i]=pntPntsTmp[i];
    }
    return pntPoints;
}


inline labelListList getPointPoints(size_t nPoints, const std::vector<facePiece>& facess)
{
    labelDynamicListList pntPntsTmp(nPoints);
    labelListList pntPoints(nPoints);
    for(auto& faces:facess)for_(faces,fI)
    {
        const face & f=faces[fI];
        for(auto pItr = f.cbegin(); pItr != f.cend(); ++pItr)
        {
           appendUnique(pntPntsTmp[ *pItr ], *nextCircIter(f,pItr));
           appendUnique(pntPntsTmp[ *pItr ], *prevCircIter(f,pItr));
        }
    }

    for_i(pntPoints)
    {
        pntPoints[i]=pntPntsTmp[i];
    }
    return pntPoints;
}

inline labelListList pointFaces(size_t nPoints, const facePiece& faces) {
    labelDynamicListList pointFacesTmp(nPoints);
    labelListList pFaces(nPoints);
    for_(faces,fI)
    {
        const face f=faces[fI];
        for(auto pI:f) appendUnique(pointFacesTmp[pI], fI);
    }

    for_i(pFaces)
    {
        pFaces[i]=pointFacesTmp[i];
    }
    return pFaces;
}
inline labelListList edgeFaces(const ints& myPPoints, const facePiece& faces, const ints& myPFaces) {
    labelDynamicListList edgFacsTmp(myPPoints.size());
    labelListList edgFacs(myPPoints.size());
    for_(myPFaces,mfI)
    {
        const face & f=faces[myPFaces[mfI]];
        for(auto pI:f)
        {
			for_(myPPoints,neIp)
				if(myPPoints[neIp]==pI)
				{ appendUnique(edgFacsTmp[neIp], myPFaces[mfI]); break; }
        }
    }

    for_i(edgFacs)
    {
        edgFacs[i]=edgFacsTmp[i];
    }
    return edgFacs;
}
inline dbl3 areax2(const face& f, const piece<point>& points) {
	dbl3 diag = points[f[2]]- points[f[0]];	return ((points[f[1]]-points[f[0]])^diag) + (diag^(points[f[3]]-points[f[0]]));
}
inline dbl3 normal(const face& f, const piece<point>& points) {	dbl3 aa = areax2(f,points);	return aa/mag(aa); }
inline dbl3 centre(const face& f, const piece<point>& points) {	return 0.25*(points[f[0]] + points[f[1]] + points[f[2]] + points[f[3]]); }
inline dbl3s faceCentres(const facePiece& facez, const piece<point>& points) {
	dbl3s res(facez.size());
	for(size_t fi=0; fi<facez.size(); ++fi) res[fi]=centre(facez[fi],points);
	return res;
}
inline vars<dbl3s> faceCentres(const std::vector<facePiece>& facezs, const piece<point>& points){
	vars<dbl3s> res(facezs.size());
	for(size_t iz=0;iz<facezs.size();++iz) res[iz]=faceCentres(facezs[iz],points);
	return res;
}
