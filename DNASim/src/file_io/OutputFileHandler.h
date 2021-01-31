// OutputFileHandler class
// Nicolas Clauvelin


// output (writing) file handler implementation
//
// copy constructor and operator are deleted because stream are not copiable 


#ifndef DNASim_OutputFileHandler_h
#define DNASim_OutputFileHandler_h


#include "file_io/FileHandler.h"


namespace DNASim {
	
	
	class EnhancedString;


	class OutputFileHandler final : public FileHandler {
		
		
	public:
		
		// constructors
		// default constructor creates a non-initialized file handler
		OutputFileHandler() = default;
		OutputFileHandler(const std::string& filename, const std::string& path)
        : FileHandler(filename, path), m_output() {};
		OutputFileHandler(const OutputFileHandler& ofh) = delete;
        OutputFileHandler(OutputFileHandler&& ofh) = default;
		~OutputFileHandler();

        // copy and move operators
        OutputFileHandler& operator=(const OutputFileHandler& ofh) = delete;
        OutputFileHandler& operator=(OutputFileHandler&& ofh) = default;
		
		// file properties accessors
		bool is_open() const override;
		bool is_eof() const override;
		bool is_good() const;
		
		// file methods
		bool open() override;
		bool open_in_append_mode();
		bool close() override;
		
		// ofstream accessor
		std::ofstream& stream();
		
		// writing method
		// these methods require the file to be opened
		bool write_line(const EnhancedString& line);
		bool write_line(const std::string& line);
		bool write_block_of_lines(const std::vector<EnhancedString>& lines);
		bool write_stream(const std::stringstream& stream);
		
		
	private:
		
		std::ofstream	m_output;
		
		
	};
	
	
}


#endif	// DNASim_OutputFileHandler_h
