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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::pimpleControl::nCorrPIMPLE() const
{
	return nCorrPIMPLE_;
}


inline Foam::label Foam::pimpleControl::nCorrPiso() const
{
	return nCorrPISO_;
}


inline Foam::label Foam::pimpleControl::corrPiso() const
{
	return corrPISO_;
}


inline bool Foam::pimpleControl::correct()
{
	corrPISO_++;

	if (debug)
	{
		Info<< algorithmName_ << " correct: corrPiso = " << corrPISO_ << endl;
	}

	if (corrPISO_ <= nCorrPISO_)
	{
		return true;
	}
	else
	{
		corrPISO_ = 0;
		return false;
	}
}


inline bool Foam::pimpleControl::storeInitialResiduals() const
{
	// Start from second PIMPLE iteration
	return (corr_ == 2) && (corrPISO_ == 0) && (corrNonOrtho_ == 0);
}


inline bool Foam::pimpleControl::firstIter() const
{
	return corr_ == 1;
}


inline bool Foam::pimpleControl::finalIter() const
{
	return converged_ || (corr_ == nCorrPIMPLE_);
}


inline bool Foam::pimpleControl::finalInnerIter() const
{
	return
	   finalIter()
	&& corrPISO_ == nCorrPISO_
	&& corrNonOrtho_ == nNonOrthCorr_ + 1;
}


inline bool Foam::pimpleControl::turbCorr() const
{
	return !turbOnFinalIterOnly_ || finalIter();
}


// ************************************************************************* //
