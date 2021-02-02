// OutputFileHandler class
// Nicolas Clauvelin


#include "file_io/EnhancedString.h"
#include "file_io/FileException.h"
#include "file_io/OutputFileHandler.h"


namespace DNASim {


	// class destructor
	OutputFileHandler::~OutputFileHandler() {
		close();
	};


	// stream opening checking method
	bool OutputFileHandler::is_open() const {
		if (m_output.is_open())
			return true;
		return false;
	};


	// stream EOF checking method
	bool OutputFileHandler::is_eof() const {
		// we always return false since output file are never EOF
		return false;
	};


	// stream condition checking method
	bool OutputFileHandler::is_good() const {
		if (m_output.good())
			return true;
		return false;
	};


	// file opening method
	bool OutputFileHandler::open() {
		
		//exceptions
		m_output.exceptions(std::ofstream::badbit | std::ofstream::failbit);
		
		// create directory if needed
		// unix only implementation
		struct stat st;
		if (stat(file_path().c_str(), &st) != 0) {
			Integer chk = mkdir(file_path().c_str(), 0755);
			if (chk != 0) {
				FileException fex("OutputFileHandler::open() "
								  "= failed to create directory for the output "
								  "file:");
				fex.add_info(full_file_name());
				fex.report(std::cout);
				exit(-1);
			};
		};
		
		// file opening
		try {
			m_output.open(full_file_name().c_str(), std::ios_base::out);
		}
		
		// exception catching
		catch (std::ofstream::failure& e) {
			FileException fex("OutputFileHandler::open() = "
							  "failed to open file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// file append mode opening method
	bool OutputFileHandler::open_in_append_mode() {
		
		//exceptions
		m_output.exceptions(std::ofstream::badbit);
		
		// file opening
		try {
			m_output.open(full_file_name().c_str(),
						  std::ios_base::out | std::ios_base::app);
		}
		
		// exception catching
		catch (std::ofstream::failure& e) {
			FileException fex("OutputFileHandler::open_in_append_mode() = "
							  "failed to open file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// file closing method
	bool OutputFileHandler::close() {
		
		// file closing
		try {
			if (is_open())
				m_output.close();
		}
		
		// exception catching
		catch (std::ofstream::failure& e) {
			FileException fex("OutputFileHandler::close() = "
							  "failed to close file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// output stream accessor
	std::ofstream& OutputFileHandler::stream() {
		
		return m_output;
		
	};


	// single line writing method
	bool OutputFileHandler::write_line(const EnhancedString& line) {
		
		// write line
		try {
			
			// check if is open
			if (!is_open()) {
				FileException fex("OutputFileHandler::write_line() = "
								  "failed to open file:");
				fex.add_info(full_file_name());
				throw fex;
			};
			
			m_output << line << "\n";
			
		}
		
		// exception catching
		catch (FileException& fex) {
			fex.add_info("OutputFileHandler::write_line() = "
						"failed to write file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// single line writing method
	bool OutputFileHandler::write_line(const std::string& line) {
		
		// write line
		try {
			
			// check if is open
			if (!is_open()) {
				FileException fex("OutputFileHandler::write_line() = "
								  "failed to open file:");
				fex.add_info(full_file_name());
				throw fex;
			};
			
			m_output << line << "\n";
			
		}
		
		// exception catching
		catch (FileException& fex) {
			fex.add_info("OutputFileHandler::write_line()"
						" = failed to write file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// line block writing method
	bool OutputFileHandler::write_block_of_lines(const
												 std::vector<EnhancedString>&
												 lines) {

		for (auto it = lines.begin(), end = lines.end(); it != end; ++it) {
			write_line(*it);
		};
		
		return true;
		
	};


	// stream writing method
	bool OutputFileHandler::write_stream(const std::stringstream& ostream) {
		
		// write stream
		try {
			
			// check if is open
			if (!is_open()) {
				FileException fex("OuputFileHandler::write_stream() = "
								  "failed to open file:");
				fex.add_info(full_file_name());
				throw fex;
			};
			
			m_output << ostream.str();
			
		}
		
		// exception catching
		catch (FileException& fex) {
			fex.add_info("OutputFileHandler::write_stream() = "
						 "failed to write file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


}
