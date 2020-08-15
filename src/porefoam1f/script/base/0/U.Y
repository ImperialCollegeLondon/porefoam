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

internalField   uniform (0 0 0);



boundaryField
{

	Bottom
	{
		type			zeroGradient;
	}
	Top
	{
		type			zeroGradient;
	}
	Grainwalls
	{
		type			fixedValue;
		value			uniform (0 0 0);
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
	Right
	{
		type			fixedValue;
		value			uniform (0 0 0);
	}
	Left
	{
		type			fixedValue;
		value			uniform (0 0 0);
	}
	
	default
	{
		type			fixedValue;
		value			uniform (0 0 0);
	}
	"VV.*"
	{
		type			slip;
		value			uniform (0 0 0);
	}
	"processor.*"
	{
		type			processor;
		value			uniform (0 0 0);
	}
}


// ************************************************************************* //
