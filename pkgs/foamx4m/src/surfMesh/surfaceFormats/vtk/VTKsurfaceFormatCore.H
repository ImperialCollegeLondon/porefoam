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
    Foam::fileFormats::VTKsurfaceFormatCore

Description
    Internal class used by the VTKsurfaceFormat

SourceFiles
    VTKsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef VTKsurfaceFormatCore_H
#define VTKsurfaceFormatCore_H

#include "Ostream.H"
#include "OFstream.H"
#include "MeshedSurface.H"
//#include "objectRegistry.H"


///////////////////////////////////////////
// this is a copy-past from SiR.h
///////////////////////////////////////////

#include <unordered_map>
#include <functional>
using deleter_t = std::function<void(void *)>;
using unique_void_ptr = std::unique_ptr<void, deleter_t>;
template<typename T>  auto deleter_silent (void const * data) -> void {  T const* p = static_cast<T const*>(data);   if(p) delete p; }
template<typename T>  auto deleter_verbose(void const * data) -> void {  T const* p = static_cast<T const*>(data);   if(p) { delete p;  std::cout<<" @"<<p<<" deleted."<<std::endl; } }
//template<typename T>  auto make_unique_void(T * ptr) -> unique_void_ptr { return unique_void_ptr(ptr, &deleter_verbose<T>); }

using objRegistry = std::unordered_map<std::string,unique_void_ptr>;
inline void*	dbget(objRegistry& stor, const std::string& var)
	{ auto vptr=stor.find(var); if(vptr==stor.end()) { if(var!="skip") cout<<"\n *** Error "<<var<<" not found *** n_stored:"<<stor.size()<<std::endl; return nullptr;} return vptr->second.get(); }
template<typename T>	T& dbset(objRegistry& stor, const std::string& var, T* ptr)
	{ (cout<<" keeping "<<var<<":").flush(); cout<<stor.emplace(std::make_pair(var, unique_void_ptr(ptr, &deleter_verbose<T>))).second; return *ptr;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


class VTKsurfaceFormatCore
{
protected:

    // Protected Static Member Functions


    // Data

        //- Cell based fields
        objRegistry cellData_;

        //- Point based fields
        objRegistry pointData_;

        //- Other fields
        objRegistry otherData_;



        //- Cell based fields
        const objRegistry& cellData() const
        {
            return cellData_;
        }

        objRegistry& cellData()
        {
            return cellData_;
        }

        //- Point based fields
        const objRegistry& pointData() const
        {
            return pointData_;
        }

        objRegistry& pointData()
        {
            return pointData_;
        }

        //- Other fields
        const objRegistry& otherData() const
        {
            return otherData_;
        }

        objRegistry& otherData()
        {
            return otherData_;
        }

	objRegistry& selectRegistry
	(
		const string& readMode
	)
	{
		Info<<" ***readMode: "<<readMode<<"*** ";
		if (readMode == "CELL_DATA")
		{
			return cellData_;
		}
		else if (readMode == "POINT_DATA")
		{
			return pointData_;
		}
		else
		{
			return otherData_;
		}
	}









        //- Write header information with points
        static void writeHeader
        (
			Ostream&,
			const pointField&
        );

        //- Write regions (zones) information as CellData
		static void writeTail(Ostream&, const UList<surfZone>&);

        //- Write regions (zones) information as CellData
		static void writeTail(Ostream&, const UList<label>& zoneIds);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
