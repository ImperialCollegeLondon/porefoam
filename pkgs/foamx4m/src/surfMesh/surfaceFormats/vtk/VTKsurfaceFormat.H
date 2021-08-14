/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fileFormats::VTKsurfaceFormat

Description
    Read/write VTK legacy format (ASCII) for surfaces.

    The output is never sorted by zone.

SourceFiles
    VTKsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKsurfaceFormat_H
#define VTKsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "VTKsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


template<class Face>
class VTKsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTKsurfaceFormatCore
{
    // Private Member Functions

        //- Write polygons
		static void writeHeaderPolygons(Ostream&, const UList<Face>&);

        //- No copy construct
        VTKsurfaceFormat(const VTKsurfaceFormat<Face>&) = delete;

        //- No copy assignment
        void operator=(const VTKsurfaceFormat<Face>&) = delete;


public:

    // Constructors

        //- Construct from file name
        VTKsurfaceFormat(const fileName& filename);


		static autoPtr<MeshedSurface<Face> > New(const fileName& name)
		{
			return autoPtr<MeshedSurface<Face> >
			(
				new VTKsurfaceFormat(name)
			);
		}

    //- Destructor
    virtual ~VTKsurfaceFormat() = default;


    // Static Member Functions

    // Member Functions

		// Write

		//- Write surface mesh components by proxy
		static void write(Ostream&, const MeshedSurfaceProxy<Face>&);

		static void write(const fileName&, const MeshedSurfaceProxy<Face>&);

		//- Write UnsortedMeshedSurface, the output remains unsorted
		static void write
		(
			const fileName&, const UnsortedMeshedSurface<Face>&
		);


		//- Write object
		virtual void write(Ostream& os) const
		{
			write(os, MeshedSurfaceProxy<Face>(*this));
		}

        //- Read from file
        virtual bool read(const fileName& filename);

		virtual void write(const fileName& name) const
		{
			write(name, MeshedSurfaceProxy<Face>(*this));
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VTKsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
