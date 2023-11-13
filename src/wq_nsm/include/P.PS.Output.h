#pragma once

#include <iostream>

#include <fstream>

#include <sstream>






namespace NSMPS
{

	
	// D E F I N E   O U T P U T   F I L E


	struct OutputFile
	{
		std::string filename;

		std::ofstream* of;

		OutputFile(std::string filename);

		~OutputFile();
	};



}


