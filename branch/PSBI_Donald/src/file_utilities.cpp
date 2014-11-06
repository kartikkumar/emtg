#include <fstream>
#include <iostream>

#include <boost/filesystem.hpp>

#include "file_utilities.h"

namespace fs = ::boost::filesystem;

using namespace std;

namespace EMTG {namespace filesystem {
	// return the filenames of all files that have the specified extension
	// in the specified directory and all subdirectories

	void get_all_files_with_extension(const fs::path& root, const string& ext, vector<fs::path>& ret)
	{  
	  if (!fs::exists(root)) return;

	  if (fs::is_directory(root))
	  {
		fs::recursive_directory_iterator it(root);
		fs::recursive_directory_iterator endit;
		while(it != endit)
		{
			fs::path file(*it);
		
			if (fs::is_regular_file(file) && file.extension() == ext)
			ret.push_back(file.filename());
		  ++it;
		}
	  }
	}

	int read_sequence_file(string filename, vector< vector<int> > sequence_database, int sequence_length)
	{
		ifstream inputfile(filename.c_str());

		int linenumber = 0;
		string peek;
		int value;
		char dump_buffer[1024];

		if (!inputfile.is_open())
			return 2;

		while (!inputfile.eof()) {
			peek = inputfile.peek();
			if (peek == "#" || peek == "\r" || peek == "\n") {
				//comment or blank line, do not parse
				inputfile.getline(dump_buffer, 1024);
				++linenumber;
			}
			else
			{
				vector<int> temp_sequence;
				for (int entry = 0; entry < sequence_length; ++entry)
				{
					inputfile >> value;
					temp_sequence.push_back(value);
				}

				sequence_database.push_back(temp_sequence);
			}

		}

		return 0;

	}
}}//close namespace