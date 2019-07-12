// FileException class
// Nicolas Clauvelin


// simple exception class for input/output file operations


#ifndef DNASim_FileException_h
#define DNASim_FileException_h


#include <DNASim_Includes.h>


namespace DNASim {


	class FileException {
	
	
	public:
		
		FileException() = default;
		FileException(const std::string& info);
		~FileException() = default;
		
		// add info methods
		void add_info(const std::string& info);
		
		// report method
		void report(std::ostream& output);
		
		
	private:
		
		std::vector<std::string> m_context;
		
		
	};
	
	
}


#endif	// DNASim_FileException_h
