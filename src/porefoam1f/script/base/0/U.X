FoamFile // -*- C++ -*-
{
	version		 2.0;
	format		  ascii;

	root			"";
	case			"";
	instance		"";
	local		   "";

	class		   volVectorField;
	object		  U;
}




dimensions	  [0 1 -1 0 0 0 0];

internalField   uniform (1e-6 0 0);



boundaryField
{

	Left
	{
		type			zeroGradient;
	}
	Right
	{
		type			zeroGradient;
	}
	Grainwalls
	{
		type			fixedValue;
		value		   uniform (0 0 0);
	}
	Front
	{
		type			fixedValue;
		value		   uniform (0 0 0);
	}
	Back
	{
		type			fixedValue;
		value		   uniform (0 0 0);
	}
	Top
	{
		type			fixedValue;
		value		   uniform (0 0 0);
	}
	Bottom
	{
		type			fixedValue;
		value		   uniform (0 0 0);
	}
	
	default
	{
		type			fixedValue;
		value		   uniform (0 0 0);
	}	
	"VV.*"
	{
		type			slip;
		value		   uniform (0 0 0);
	}
	"processor.*"
	{
		type			processor;
		value		   uniform (0 0 0);
	}
}


// ************************************************************************* //
