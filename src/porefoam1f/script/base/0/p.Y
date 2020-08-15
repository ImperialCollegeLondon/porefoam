FoamFile // -*- C++ -*-
{
	version	 2.0;
	format	  ascii;
	class	   volScalarField;
	location	"0";
	object	  p;
}


dimensions	  [1 -1 -2 0 0 0 0];

//dimensions	  [0 2 -2 0 0 0 0];

internalField   uniform 0;




boundaryField
{

	Grainwalls
	{
		type			zeroGradient;
	}
	Bottom
	{
		type			fixedValue;
		value uniform 1;
	}
	Top
	{
		type			fixedValue;
		value 			uniform 0;
	}
	Front
	{
		type			zeroGradient;
	}
	Back
	{
		type			zeroGradient;
	}
	Right
	{
		type			zeroGradient;
	}
	Left
	{
		type			zeroGradient;
	}
	"VV.*"
	{
		type			zeroGradient;
	}
	"processor.*"
	{
		type			processor;
		value		   uniform  0.0;
	}
}




