// FileException class
// Nicolas Clauvelin


#include "file_io/FileException.h"


namespace DNASim {


	// class constructor with initialized context
	FileException::FileException(const std::string& info) : m_context() {

		m_context.push_back(info);
		
	};


	// add info methods
	void FileException::add_info(const std::string& info) {
		
		m_context.push_back(info);
		
	};


	// report method
	void FileException::report(std::ostream& output) {
		
		output << "\n" << "[FileException]:" << "\n";
		for (Size i=0; i<m_context.size(); ++i)
			output << m_context[i] << "\n";
		output << "\n";
		
	};

	
}
