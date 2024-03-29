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

#include "fvSchemes.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::debug::debugSwitch
Foam::fvSchemes::debug
(
	"fvSchemes",
	false
);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fvSchemes::clear()
{
	ddtSchemes_.clear();
	defaultDdtScheme_.clear();
	d2dt2Schemes_.clear();
	defaultD2dt2Scheme_.clear();
	interpolationSchemes_.clear();
	defaultInterpolationScheme_.clear();
	divSchemes_.clear(); // optional
	defaultDivScheme_.clear();
	gradSchemes_.clear(); // optional
	defaultGradScheme_.clear();
	snGradSchemes_.clear();
	defaultSnGradScheme_.clear();
	laplacianSchemes_.clear(); // optional
	defaultLaplacianScheme_.clear();
	fluxRequired_.clear();
	defaultFluxRequired_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvSchemes::fvSchemes(const objectRegistry& obr)
:
	IOdictionary
	(
		IOobject
		(
			"fvSchemes",
			obr.time().system(),
			obr,
//			 IOobject::MUST_READ,
			IOobject::READ_IF_PRESENT,  // Allow default dictionary creation
			IOobject::NO_WRITE
		)
	),
	ddtSchemes_
	(
		ITstream
		(
			objectPath() + "::ddtSchemes",
			tokenList()
		)()
	),
	defaultDdtScheme_
	(
		ddtSchemes_.name() + "::default",
		tokenList()
	),
	d2dt2Schemes_
	(
		ITstream
		(
			objectPath() + "::d2dt2Schemes",
			tokenList()
		)()
	),
	defaultD2dt2Scheme_
	(
		d2dt2Schemes_.name() + "::default",
		tokenList()
	),
	interpolationSchemes_
	(
		ITstream
		(
			objectPath() + "::interpolationSchemes",
			tokenList()
		)()
	),
	defaultInterpolationScheme_
	(
		interpolationSchemes_.name() + "::default",
		tokenList()
	),
	divSchemes_
	(
		ITstream
		(
			objectPath() + "::divSchemes",
			tokenList()
		)()
	),
	defaultDivScheme_
	(
		divSchemes_.name() + "::default",
		tokenList()
	),
	gradSchemes_
	(
		ITstream
		(
			objectPath() + "::gradSchemes",
			tokenList()
		)()
	),
	defaultGradScheme_
	(
		gradSchemes_.name() + "::default",
		tokenList()
	),
	snGradSchemes_
	(
		ITstream
		(
			objectPath() + "::snGradSchemes",
			tokenList()
		)()
	),
	defaultSnGradScheme_
	(
		snGradSchemes_.name() + "::default",
		tokenList()
	),
	laplacianSchemes_
	(
		ITstream
		(
			objectPath() + "::laplacianSchemes",
			tokenList()
		)()
	),
	defaultLaplacianScheme_
	(
		laplacianSchemes_.name() + "::default",
		tokenList()
	),
	fluxRequired_
	(
		ITstream
		(
			objectPath() + "::fluxRequired",
			tokenList()
		)()
	),
	defaultFluxRequired_(false)
{
	if (!headerOk())
	{
		if (debug)
		{
			InfoIn
			(
				"fvSchemes::fvSchemes(const objectRegistry& obr)"
			)   << "fvSchemes dictionary not found.  Creating default."
				<< endl;
		}

		regIOobject::write();
	}

	read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvSchemes::read()
{
	if (regIOobject::read())
	{
		const dictionary& dict = schemesDict();

		// Persistent settings across reads is incorrect
		clear();

		if (dict.found("ddtSchemes"))
		{
			ddtSchemes_ = dict.subDict("ddtSchemes");
		}
		else if (dict.found("timeScheme"))
		{
			// For backward compatibility.
			// The timeScheme will be deprecated with warning or removed
			WarningIn("fvSchemes::read()")
				<< "Using deprecated 'timeScheme' instead of 'ddtSchemes'"
				<< nl << endl;

			word schemeName(dict.lookup("timeScheme"));

			if (schemeName == "EulerImplicit")
			{
				schemeName = "Euler";
			}
			else if (schemeName == "BackwardDifferencing")
			{
				schemeName = "backward";
			}
			else if (schemeName == "SteadyState")
			{
				schemeName = "steadyState";
			}
			else
			{
				FatalIOErrorIn("fvSchemes::read()", dict.lookup("timeScheme"))
					<< "\n    Only EulerImplicit, BackwardDifferencing and "
					   "SteadyState\n    are supported by the old timeScheme "
					   "specification.\n    Please use ddtSchemes instead."
					<< exit(FatalIOError);
			}

			ddtSchemes_.set("default", schemeName);

			ddtSchemes_.lookup("default")[0].lineNumber() =
				dict.lookup("timeScheme").lineNumber();
		}
		else
		{
			ddtSchemes_.set("default", "none");
		}

		if
		(
			ddtSchemes_.found("default")
		 && word(ddtSchemes_.lookup("default")) != "none"
		)
		{
			defaultDdtScheme_ = ddtSchemes_.lookup("default");
		}


		if (dict.found("d2dt2Schemes"))
		{
			d2dt2Schemes_ = dict.subDict("d2dt2Schemes");
		}
		else if (dict.found("timeScheme"))
		{
			// For backward compatibility.
			// The timeScheme will be deprecated with warning or removed
			WarningIn("fvSchemes::read()")
				<< "Using deprecated 'timeScheme' instead of 'd2dt2Schemes'"
				<< nl << endl;

			word schemeName(dict.lookup("timeScheme"));

			if (schemeName == "EulerImplicit")
			{
				schemeName = "Euler";
			}
			else if (schemeName == "SteadyState")
			{
				schemeName = "steadyState";
			}

			d2dt2Schemes_.set("default", schemeName);

			d2dt2Schemes_.lookup("default")[0].lineNumber() =
				dict.lookup("timeScheme").lineNumber();
		}
		else
		{
			d2dt2Schemes_.set("default", "none");
		}

		if
		(
			d2dt2Schemes_.found("default")
		 && word(d2dt2Schemes_.lookup("default")) != "none"
		)
		{
			defaultD2dt2Scheme_ = d2dt2Schemes_.lookup("default");
		}


		if (dict.found("interpolationSchemes"))
		{
			interpolationSchemes_ = dict.subDict("interpolationSchemes");
		}
		else if (!interpolationSchemes_.found("default"))
		{
			interpolationSchemes_.add("default", "linear");
		}

		if
		(
			interpolationSchemes_.found("default")
		 && word(interpolationSchemes_.lookup("default")) != "none"
		)
		{
			defaultInterpolationScheme_ =
				interpolationSchemes_.lookup("default");
		}


		if (dict.found("divSchemes"))
		{
			divSchemes_ = dict.subDict("divSchemes");
		}

		if
		(
			divSchemes_.found("default")
		 && word(divSchemes_.lookup("default")) != "none"
		)
		{
			defaultDivScheme_ = divSchemes_.lookup("default");
		}

		if (dict.found("gradSchemes"))
		{
			gradSchemes_ = dict.subDict("gradSchemes");
		}

		if
		(
			gradSchemes_.found("default")
		 && word(gradSchemes_.lookup("default")) != "none"
		)
		{
			defaultGradScheme_ = gradSchemes_.lookup("default");
		}


		if (dict.found("snGradSchemes"))
		{
			snGradSchemes_ = dict.subDict("snGradSchemes");
		}
		else if (!snGradSchemes_.found("default"))
		{
			snGradSchemes_.add("default", "corrected");
		}

		if
		(
			snGradSchemes_.found("default")
		 && word(snGradSchemes_.lookup("default")) != "none"
		)
		{
			defaultSnGradScheme_ = snGradSchemes_.lookup("default");
		}


		if (dict.found("laplacianSchemes"))
		{
			laplacianSchemes_ = dict.subDict("laplacianSchemes");
		}

		if
		(
			laplacianSchemes_.found("default")
		 && word(laplacianSchemes_.lookup("default")) != "none"
		)
		{
			defaultLaplacianScheme_ = laplacianSchemes_.lookup("default");
		}


		if (dict.found("fluxRequired"))
		{
		     fluxRequired_ = dict.subDict("fluxRequired");

			if
			(
				fluxRequired_.found("default")
			 && word(fluxRequired_.lookup("default")) != "none"
			)
			{
				defaultFluxRequired_ = Switch(fluxRequired_.lookup("default"));
			}
		}

		return true;
	}
	else
	{
		return false;
	}
}


const Foam::dictionary& Foam::fvSchemes::schemesDict() const
{
	if (found("select"))
	{
		return subDict(word(lookup("select")));
	}
	else
	{
		return *this;
	}
}


Foam::ITstream& Foam::fvSchemes::ddtScheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup ddtScheme for " << name << endl;
	}

	if (ddtSchemes_.found(name) || defaultDdtScheme_.empty())
	{
		return ddtSchemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultDdtScheme_).rewind();
		return const_cast<ITstream&>(defaultDdtScheme_);
	}
}


Foam::ITstream& Foam::fvSchemes::d2dt2Scheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup d2dt2Scheme for " << name << endl;
	}

	if (d2dt2Schemes_.found(name) || defaultD2dt2Scheme_.empty())
	{
		return d2dt2Schemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultD2dt2Scheme_).rewind();
		return const_cast<ITstream&>(defaultD2dt2Scheme_);
	}
}


Foam::ITstream& Foam::fvSchemes::interpolationScheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup interpolationScheme for " << name << endl;
	}

	if
	(
		interpolationSchemes_.found(name)
	 || defaultInterpolationScheme_.empty()
	)
	{
		return interpolationSchemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultInterpolationScheme_).rewind();
		return const_cast<ITstream&>(defaultInterpolationScheme_);
	}
}


Foam::ITstream& Foam::fvSchemes::divScheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup divScheme for " << name << endl;
	}

	if (divSchemes_.found(name) || defaultDivScheme_.empty())
	{
		return divSchemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultDivScheme_).rewind();
		return const_cast<ITstream&>(defaultDivScheme_);
	}
}


Foam::ITstream& Foam::fvSchemes::gradScheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup gradScheme for " << name << endl;
	}

	if (gradSchemes_.found(name) || defaultGradScheme_.empty())
	{
		return gradSchemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultGradScheme_).rewind();
		return const_cast<ITstream&>(defaultGradScheme_);
	}
}


Foam::ITstream& Foam::fvSchemes::snGradScheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup snGradScheme for " << name << endl;
	}

	if (snGradSchemes_.found(name) || defaultSnGradScheme_.empty())
	{
		return snGradSchemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultSnGradScheme_).rewind();
		return const_cast<ITstream&>(defaultSnGradScheme_);
	}
}


Foam::ITstream& Foam::fvSchemes::laplacianScheme(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup laplacianScheme for " << name << endl;
	}

	if (laplacianSchemes_.found(name) || defaultLaplacianScheme_.empty())
	{
		return laplacianSchemes_.lookup(name);
	}
	else
	{
		const_cast<ITstream&>(defaultLaplacianScheme_).rewind();
		return const_cast<ITstream&>(defaultLaplacianScheme_);
	}
}


void Foam::fvSchemes::setFluxRequired(const word& name) const
{
	if (debug)
	{
		Info<< "Setting fluxRequired for " << name << endl;
	}

	fluxRequired_.add(name, true, true);
}


bool Foam::fvSchemes::fluxRequired(const word& name) const
{
	if (debug)
	{
		Info<< "Lookup fluxRequired for " << name << endl;
	}

	if (fluxRequired_.found(name))
	{
		return true;
	}
	else
	{
		return defaultFluxRequired_;
	}
}

bool Foam::fvSchemes::writeData(Ostream& os) const
{
	// Write dictionaries
	os << nl << "ddtSchemes";
	ddtSchemes_.write(os, true);

	os << nl << "d2dt2Schemes";
	d2dt2Schemes_.write(os, true);

	os << nl << "interpolationSchemes";
	interpolationSchemes_.write(os, true);

	os << nl << "divSchemes";
	divSchemes_.write(os, true);

	os << nl << "gradSchemes";
	gradSchemes_.write(os, true);

	os << nl << "snGradSchemes";
	snGradSchemes_.write(os, true);

	os << nl << "laplacianSchemes";
	laplacianSchemes_.write(os, true);

	os << nl << "fluxRequired";
	fluxRequired_.write(os, true);

	return true;
}


// ************************************************************************* //
