// FileHandler class
// Nicolas Clauvelin


// abstract class for input/output file handlers


#ifndef DNASim_FileHandler_h
#define DNASim_FileHandler_h


#include <DNASim_Includes.h>


namespace DNASim {
	
	
	class FileException;
	class EnhancedString;


	class FileHandler {
		
		
	public:
		
		// virtual destructor
		virtual ~FileHandler();
		
		// file handler properties accessor
		bool is_empty() const;
		virtual bool is_open() const =0;
		virtual bool is_eof() const =0;
		
		// file path methods
		const std::string current_directory() const;
		const std::string& file_path() const;
		void set_file_path(const std::string& path);
		
		// file methods
		virtual bool open() =0;
		virtual bool close() =0;
		const std::string& file_name() const;
		const std::string full_file_name() const;
		void set_file_name(const std::string& filename);
		
		
	protected:
		
		// constructors
		// default constructor creates an empty handler
		// if not initialized the file path is set to DOT_DIRECTORY
		FileHandler() = default;
		FileHandler(const std::string& filename);	// file name initialization
		FileHandler(const std::string& filename,	// file name and path init
					const std::string& path);
		FileHandler(const FileHandler& fh) = default;
        FileHandler(FileHandler&& fh) = default;

        // copy and move operators
        FileHandler& operator=(const FileHandler& fh) = default;
        FileHandler& operator=(FileHandler&& fh) = default;
		
		
	private:
		
		std::string	m_file_name;
		std::string	m_file_path;
		
		
	};
	
	
}


#endif	// DNASim_FileHandler_h
