//filesystem utilities for EMTG

#ifndef _EMTG_FILESYSTEM
#define _EMTG_FILESYSTEM

#define BOOST_FILESYSTEM_VERSION 3
#define BOOST_FILESYSTEM_NO_DEPRECATED 

#include <vector>
#include <string>

#include <boost/filesystem.hpp>


namespace fs = ::boost::filesystem;
using namespace std;

namespace EMTG { namespace filesystem {

	void get_all_files_with_extension(const fs::path& root, const string& ext, vector<fs::path>& ret);
	int read_sequence_file(string filename, vector< vector<int> > sequence_database, int sequence_length);
}}

#endif //_EMTG_FILESYSTEM