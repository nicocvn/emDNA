// FileHandler class
// Nicolas Clauvelin


#include <FileException.h>
#include <EnhancedString.h>
#include <FileHandler.h>


// file i/o constants
#define	DOT_DIRECTORY '.'
#define	PATH_DELIMITER '/'


namespace DNASim {


	// class constructor with initialized filename
	FileHandler::FileHandler(const std::string& filename) :
	m_file_name(), m_file_path() {
		
		set_file_name(filename);
		
		// assume the file to be in the current directory
		set_file_path(std::string((const char* const)DOT_DIRECTORY));
		
	};


	// class constructor with initialized file name and path
	FileHandler::FileHandler(const std::string& filename,
							 const std::string& path) :
	m_file_name(filename), m_file_path(path) {};


	// return true if both file name and path are empty string
	bool FileHandler::is_empty() const {
		
		if (m_file_name.empty() == 0 &&
			m_file_path.empty() == 0)
			return true;
		
		return false;
		
	};


    // class destructor
    FileHandler::~FileHandler() = default;


	// return the current working directory
	const std::string FileHandler::current_directory() const {
#define BUFFER_SIZE 2048
		char cwd[BUFFER_SIZE];
		char* ret = getcwd(cwd, BUFFER_SIZE);
        assert(ret != nullptr);
#undef BUFFER_SIZE
		return std::string(cwd);
	};


	// file path accessor
	const std::string& FileHandler::file_path() const {
		return m_file_path;
	};


	// file path modifier
	void FileHandler::set_file_path(const std::string& path) {
		m_file_path = path;
	};


	// file name accessor
	const std::string& FileHandler::file_name() const {
		return m_file_name;
	};


	// full file name accessor
	const std::string FileHandler::full_file_name() const {
		
		// string for storing the full path name
		std::string fullpath;
		
		// path last caracter extraction
		std::string::const_iterator StringIt = m_file_path.end();
		--StringIt;
		const char last = *StringIt;
		
		// full path building; first add a /
        const char delm = PATH_DELIMITER;
		if (strcmp(&last, &delm) == 0)
			fullpath = m_file_path;
		else
			fullpath = m_file_path + "/";
		fullpath += m_file_name;
		
		return fullpath;
		
	};


	// file name modifier
	void FileHandler::set_file_name(const std::string& filename) {
		m_file_name = filename;
	};


}


#undef DOT_DIRECTORY
#undef PATH_DELIMITER
