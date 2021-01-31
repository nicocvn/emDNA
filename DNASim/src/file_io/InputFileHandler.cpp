// InputFileHandler class
// Nicolas Clauvelin


#include "file_io/EnhancedString.h"
#include "file_io/FileException.h"
#include "file_io/InputFileHandler.h"


namespace DNASim {


	// class destructor
	InputFileHandler::~InputFileHandler() {
		close();
	};


	// stream opening checking method
	bool InputFileHandler::is_open() const {
		if (m_input.is_open())
			return true;
		return false;
	};


	// stream EOF checking method
	bool InputFileHandler::is_eof() const {
		if (m_input.eof())
			return true;
		return false;
	};


	// stream condition checking method
	bool InputFileHandler::is_good() const {
		if (m_input.good())
			return true;
		return false;
	};


	// file opening method
	bool InputFileHandler::open() {
		
		// exceptions
		// failbit is not set in order to detect EOF
		m_input.exceptions(std::ifstream::badbit);
		
		// file opening
		try {
			m_input.open(full_file_name().c_str(), std::ios_base::in);
		}
		
		// exception catching
		catch (std::ifstream::failure& e) {
			FileException fex("InputFileHandler::open() = "
							  "failed to open file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// file closing method
	bool InputFileHandler::close() {
		
		// file closing
		try {
			if (is_open())
				m_input.close();
		}
		
		// exception catching
		catch (std::ifstream::failure& e) {
			FileException fex("InputFileHandler::close() = "
							  "failed to close file:");
			fex.add_info(full_file_name());
			fex.report(std::cout);
			exit(-1);
		};
		
		return true;
		
	};


	// read line method
	bool InputFileHandler::read_line(EnhancedString& line) {
		
		// read line
		std::string sline;
		try {
			
			// check if open
			if (!is_open()) {
				FileException fex("InputFileHandler::read_line() = "
								  "failed to open file:");
				fex.add_info(full_file_name());
				fex.report(std::cout);
				exit(-1);
			};
			
			getline(m_input, sline);
			if (is_eof())
				return false;
			
		}
		
		// exception catching
		catch (std::ifstream::failure& e) {
			FileException fex("InputFileHandler::read_line() = "
							  "failed to read file:");
			fex.add_info(full_file_name());
			throw fex;
		};
		
		// convert to EnhancedString
		line = sline;
		
		return true;
		
	};


	// file reading method
	void InputFileHandler::read_file(std::vector<EnhancedString>& filelines) {
		
		try {
		
			// file opening
			open();
		
			// clear container
			filelines.clear();
		
			// reading loop
			EnhancedString line;
			while(read_line(line))
				filelines.push_back(line);
		
			// file closing
			close();
			
		}
		
		// exception
		catch (FileException& fex) {
			fex.add_info("InputFileHandler::read_file() = failed to read file");
			fex.report(std::cout);
			exit(-1);
		};
		
	};


	// single key values file reading method
	void
	InputFileHandler::read_single_key_value_file(std::map<std::string,
												 EnhancedString>& keys,
												 const char& separator,
												 const char& commentchar) {
		
		// clear container
		keys.clear();
		
		// regular file reading
		std::vector<EnhancedString> lines;
		try { read_file(lines); }
		catch (FileException& fex) {
			fex.add_info("InputFileHandler::read_single_key_values_file()");
			fex.report(std::cout);
			exit(-1);
		};
		
		// parsing
		for (Size i=0; i<lines.size(); ++i) {
			
			// reject comment lines
			// strcmp returns 0 if strings are different
			const char ch = (std::string(lines[i]).data())[0];
			if (ch != commentchar) {
				
				// tokenizing
				std::vector<EnhancedString> tks;
				bool tk = lines[i].tokenize(tks, separator);
			
				// insert key with value if correctly formated
				// (conversion to std::string is implicit due to cast operator)
				if (tk && tks.size()==2)
					keys.
                    insert(std::
                           pair<std::string,EnhancedString>(std::string(tks[0]),
                                                            tks[1]));
			
			};
			
		};
		
	};


    // file reading static method
    std::vector<std::string>
    InputFileHandler::read_file(const std::string& fname,
                                const std::string& dir) {

        // handler
        InputFileHandler ifh(fname, dir);

        // reading
        std::vector<EnhancedString> es_lines;
        ifh.read_file(es_lines);

        // conversion
        std::vector<std::string> lines;
        lines.reserve(es_lines.size());
        for (auto it = es_lines.begin(), end = es_lines.end(); it != end; ++it)
            lines.push_back(std::string(*it));

        return lines;

    };


    // key/value file reading static method
    std::map<std::string,std::string>
    InputFileHandler::read_key_value_file(const std::string& fname,
                                          const std::string& dir,
                                          const char& sep,
                                          const char& comment) {

        // handler
        InputFileHandler ifh(fname, dir);

        // reading
        std::map<std::string, EnhancedString> keys;
        ifh.read_single_key_value_file(keys, sep, comment);

        // conversion
        std::map<std::string,std::string> string_keys;
        for (auto it = keys.begin(), end = keys.end(); it != end; ++it) {
            std::pair<std::string,std::string>
            p(EnhancedString::remove_spaces(it->first),
              EnhancedString::remove_spaces(std::string(it->second)));
            string_keys.insert(p);
        };

        return string_keys;

    };


}
