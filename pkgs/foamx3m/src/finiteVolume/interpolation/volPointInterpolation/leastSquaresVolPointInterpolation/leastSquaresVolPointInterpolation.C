/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(leastSquaresVolPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  void leastSquaresVolPointInterpolation::calcA(List<scalarSquareMatrix>& A) const
  {
	//Info << "leastSquaresVolPointInterpolation calcA" << endl;

	const fvMesh& mesh = mesh_;
	const pointField& points = mesh.points();

	//- construct 4x4 A matrix for each point
	//List<scalarSquareMatrix>& A = A_;

	//- populate A matrix
	forAll(points, pointI)
	  {
	const labelList& pointCells = mesh.pointCells()[pointI];

	//- this component of matrix does not depend on coordinates
	A[pointI][3][3] = pointCells.size();

	//- fill the A matrices
	forAll(pointCells, pointCelli)
	  {
		const label& celli = pointCells[pointCelli];

		const scalar& x = mesh.C()[celli].component(vector::X);
		const scalar& y = mesh.C()[celli].component(vector::Y);
		const scalar& z = mesh.C()[celli].component(vector::Z);

		A[pointI][0][0] += x*x;
		A[pointI][0][1] += x*y;
		A[pointI][0][2] += x*z;
		A[pointI][0][3] += x;

		A[pointI][1][0] += x*y;
		A[pointI][1][1] += y*y;
		A[pointI][1][2] += y*z;
		A[pointI][1][3] += y;

		A[pointI][2][0] += x*z;
		A[pointI][2][1] += y*z;
		A[pointI][2][2] += z*z;
		A[pointI][2][3] += z;

		A[pointI][3][0] += x;
		A[pointI][3][1] += y;
		A[pointI][3][2] += z;
		//A[pointI][3][3] = pointCells.size(); // set above
	  }
	  }

	//- for boundary points we will include the surrounding face centres
	forAll(mesh.boundary(), patchi)
	  {
	const vectorField& faceCentres = mesh.boundaryMesh()[patchi].faceCentres();
	const labelListList& pointFaces = mesh.boundaryMesh()[patchi].pointFaces();

	if(mesh.boundary()[patchi].coupled()) //- for proc boundaries
	  {
		//- for coupled patches we will use the values at the neighbourField cell centres and we will
		//- not use the boundary face values
		//- neighbour cell centre are equal to the faceCell centres plus the delta vector
		vectorField pDelta = mesh.boundary()[patchi].delta();
		vectorField faceCellC(faceCentres.size(), vector::zero);
		forAll(faceCentres, faceI)
		  {
		label celli = mesh.boundaryMesh()[patchi].faceCells()[faceI];
		faceCellC[faceI] = mesh.C()[celli];
		  }
		vectorField neiCellC = faceCellC + pDelta;

		forAll(pointFaces, pointI)
		  {
		forAll(pointFaces[pointI], pointFacei)
		  {
			label neiCelli = pointFaces[pointI][pointFacei];
			const scalar& x = neiCellC[neiCelli].component(vector::X);
			const scalar& y = neiCellC[neiCelli].component(vector::Y);
			const scalar& z = neiCellC[neiCelli].component(vector::Z);

			label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointI];

			A[globalPointi][0][0] += x*x;
			A[globalPointi][0][1] += x*y;
			A[globalPointi][0][2] += x*z;
			A[globalPointi][0][3] += x;

			A[globalPointi][1][0] += x*y;
			A[globalPointi][1][1] += y*y;
			A[globalPointi][1][2] += y*z;
			A[globalPointi][1][3] += y;

			A[globalPointi][2][0] += x*z;
			A[globalPointi][2][1] += y*z;
			A[globalPointi][2][2] += z*z;
			A[globalPointi][2][3] += z;

			A[globalPointi][3][0] += x;
			A[globalPointi][3][1] += y;
			A[globalPointi][3][2] += z;
			A[globalPointi][3][3] += 1; // = pointCells.size();
		  }
		  }
	  }
	else
	  {
		//- each point must use at least 4 neighbouring locations otherwise A is singular
		//- and simpleMatrix will cannot invert it
		//- therefore empty patches values are included to make sure A is not singular
		forAll(pointFaces, pointI)
		  {
		label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointI];

		forAll(pointFaces[pointI], pointFacei)
		  {
			//- fix: use pointFace not face philipc
			label faceI = pointFaces[pointI][pointFacei];
			const scalar& x = faceCentres[faceI].component(vector::X);
			const scalar& y = faceCentres[faceI].component(vector::Y);
			const scalar& z = faceCentres[faceI].component(vector::Z);

			A[globalPointi][0][0] += x*x;
			A[globalPointi][0][1] += x*y;
			A[globalPointi][0][2] += x*z;
			A[globalPointi][0][3] += x;

			A[globalPointi][1][0] += x*y;
			A[globalPointi][1][1] += y*y;
			A[globalPointi][1][2] += y*z;
			A[globalPointi][1][3] += y;

			A[globalPointi][2][0] += x*z;
			A[globalPointi][2][1] += y*z;
			A[globalPointi][2][2] += z*z;
			A[globalPointi][2][3] += z;

			A[globalPointi][3][0] += x;
			A[globalPointi][3][1] += y;
			A[globalPointi][3][2] += z;
			A[globalPointi][3][3] += 1; // = pointCells.size();
		  }
		  }
	  } //- end of else
	  } //- end of forAll boundary
  }


  void leastSquaresVolPointInterpolation::calcB(List<Field<vector> >& B, const GeometricField<vector, fvPatchField, volMesh>& vf) const
  {
	//Info << "leastSquaresVolPointInterpolation calcB" << endl;

	const fvMesh& mesh = mesh_;
	const pointField& points = mesh.points();

	for (direction compi = 0; compi < 3; compi++)
	  {
	forAll(points, pointI)
	  {
		const labelList& pointCells = mesh.pointCells()[pointI];
		//Info << "\npointCells " << pointCells << endl;
		forAll(pointCells, pointCelli)
		  {
		const label& celli = pointCells[pointCelli];
		//Info << "celli " << celli << ", C is " << mesh.C()[celli] << ", mesh.C().size() " << mesh.C().size() << endl;
		//Info << "mesh.C() is " << mesh.C().internalField() << endl;
		const scalar& x = mesh.C()[celli].component(vector::X);
		const scalar& y = mesh.C()[celli].component(vector::Y);
		const scalar& z = mesh.C()[celli].component(vector::Z);

		const scalar& phiCompi = vf.internalField()[celli].component(compi);

		B[pointI][0].component(compi) += phiCompi*x;
		B[pointI][1].component(compi) += phiCompi*y;
		B[pointI][2].component(compi) += phiCompi*z;
		B[pointI][3].component(compi) += phiCompi;
		  }
	  }

	//- for boundary points we will include the surrounding face centres
	forAll(mesh.boundary(), patchi)
	  {
		const vectorField& faceCentres = mesh.boundaryMesh()[patchi].faceCentres();
		const labelListList& pointFaces = mesh.boundaryMesh()[patchi].pointFaces();
		const labelList& faceCells = mesh.boundaryMesh()[patchi].faceCells();

		//- fix: do not calculate B for empty patches - philipc
		if(mesh.boundary()[patchi].coupled())
		  {
		//- for coupled patches we will use the values at the neighbourField cell centres and we will
		//- not use the boundary face values
		//- neighbour cell centre are equal to the faceCell centres plus the delta vector
		vectorField pDelta = mesh.boundary()[patchi].delta();
		vectorField faceCellC(faceCentres.size(), vector::zero);
		forAll(faceCentres, faceI)
		  {
			label celli = mesh.boundaryMesh()[patchi].faceCells()[faceI];
			faceCellC[faceI] = mesh.C()[celli];
		  }
		vectorField neiCellC = faceCellC + pDelta;

		vectorField phiNeiField = vf.boundaryField()[patchi].patchNeighbourField();

		forAll(pointFaces, pointI)
		  {
			forAll(pointFaces[pointI], pointFacei)
			  {
			label neiCelli = pointFaces[pointI][pointFacei];
			const scalar& x = neiCellC[neiCelli].component(vector::X);
			const scalar& y = neiCellC[neiCelli].component(vector::Y);
			const scalar& z = neiCellC[neiCelli].component(vector::Z);

			label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointI];

			//- this is the value of phi at the cell centre in the neighbour (i.e. across the interface)
			scalar phiCompi = phiNeiField[neiCelli].component(compi);

			B[globalPointi][0].component(compi) += phiCompi*x;
			B[globalPointi][1].component(compi) += phiCompi*y;
			B[globalPointi][2].component(compi) += phiCompi*z;
			B[globalPointi][3].component(compi) += phiCompi;
			  }
		  }
		  }
		else
		  {
		//- each point must use at least 4 neighbouring locations otherwise A is singular
		//- and simpleMatrix will cannot invert it
		//- therefore empty patches values are included to make sure A is not singular
		forAll(pointFaces, pointI)
		  {
			forAll(pointFaces[pointI], pointFacei)
			  {
			//- fix: use pointFace not face philipc
			label faceI = pointFaces[pointI][pointFacei];
			const scalar& x = faceCentres[faceI].component(vector::X);
			const scalar& y = faceCentres[faceI].component(vector::Y);
			const scalar& z = faceCentres[faceI].component(vector::Z);

			label globalPointi = mesh.boundaryMesh()[patchi].meshPoints()[pointI];

			scalar phiCompi = 0.0;
			if(mesh.boundary()[patchi].type() == "empty")
			  {
				//- use faceCell value for empty because empty patches do not store any values
				const label& ci = faceCells[faceI];
				phiCompi = vf.internalField()[ci].component(compi);
			  }
			else
			  {
				phiCompi = vf.boundaryField()[patchi][faceI].component(compi);
			  }

				B[globalPointi][0].component(compi) += phiCompi*x;
			B[globalPointi][1].component(compi) += phiCompi*y;
			B[globalPointi][2].component(compi) += phiCompi*z;
			B[globalPointi][3].component(compi) += phiCompi;
			  }
		  }
		  }
	  } //- end of forAll boundary
	  } //- end of for all components
  }


  void leastSquaresVolPointInterpolation::interpolate
  (
   const GeometricField<vector, fvPatchField, volMesh>& vf,
   GeometricField<vector, pointPatchField, pointMesh>& pf //Field<vector>& pf
   ) const
  {
	//Info << "Interpolating cell to point using leastSquaresVolPointInterpolation" << endl;

	const fvMesh& mesh = mesh_;
	const pointField& points = mesh.points();

	//- first check that point field is the correct size
	if(pf.size() != points.size())
	  {
	FatalError << "pointfield should be equal to the number of points in the fvMesh"
		   << abort(FatalError);
	  }

	//- calculate A and B vector
	List<scalarSquareMatrix> A(mesh.points().size(), scalarSquareMatrix(4, 0.0));
	calcA(A);

	List<Field<vector> > B(mesh.points().size(), Field<vector>(4, pTraits<vector>::zero));
	calcB(B, vf);

	//- solve equations for each component of each point
	forAll(points, pointI)
	  {
	Field<vector>& source = B[pointI];
	simpleMatrix<vector> leastSquaresMatrix(A[pointI], source);
	//Info << "solving equation for point " << pointI << endl;
	//Info << "A[pointI] is " << A[pointI] << ", source is " << source << endl;
	//- solve using Gauss elimination or LU decomposition with pivoting
	//Field<vector> leastSquaresSol = leastSquaresMatrix.solve();
	Field<vector> leastSquaresSol = leastSquaresMatrix.LUsolve();

	const scalar& x = mesh.points()[pointI].component(vector::X);
	const scalar& y = mesh.points()[pointI].component(vector::Y);
	const scalar& z = mesh.points()[pointI].component(vector::Z);

	//- calculate phi at vertex
	for (direction compi = 0; compi < 3; compi++)
	  {
		const scalar& a = leastSquaresSol[0].component(compi);
		const scalar& b = leastSquaresSol[1].component(compi);
		const scalar& c = leastSquaresSol[2].component(compi);
		const scalar& d = leastSquaresSol[3].component(compi);

		pf[pointI].component(compi) = a*x + b*y + c*z + d;
	  }
	  }

	//- proc patches are synchronised
	pf.correctBoundaryConditions();
  }


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

  leastSquaresVolPointInterpolation::leastSquaresVolPointInterpolation(const fvMesh& vm)
  :
	MeshObject<fvMesh, leastSquaresVolPointInterpolation>(vm),
	mesh_(vm) //,
	//A_(vm.points().size(), scalarSquareMatrix(4, 0.0)),
	//B_(vm.points().size(), Field<vector>(4, vector::zero))
  {
	//calcA();
  }

// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

  leastSquaresVolPointInterpolation::~leastSquaresVolPointInterpolation()
  {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
