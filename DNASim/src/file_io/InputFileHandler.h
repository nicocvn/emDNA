// InputFileHandler class
// Nicolas Clauvelin


// input (reading) file handler implementation
//
// copy constructor and operator are deleted because stream are not copiable


#ifndef DNASim_InputFileHandler_h
#define DNASim_InputFileHandler_h


#include <FileHandler.h>
#include <EnhancedString.h>


namespace DNASim {
	
	
	class EnhancedString;


	class InputFileHandler final : public FileHandler {
		
		
	public:
		
		// constructors
		// default constructor creates a non-initialized file handler
		InputFileHandler() = default;
		InputFileHandler(const std::string& filename, const std::string& path) :
        FileHandler(filename, path), m_input() {};
		InputFileHandler(const InputFileHandler& ifh) = delete;
        InputFileHandler(InputFileHandler&& ifh) = default;
		~InputFileHandler();

        // copy and move operators
        InputFileHandler& operator=(const InputFileHandler& ifh) = delete;
        InputFileHandler& operator=(InputFileHandler&& ifh) = default;
		
		// file properties accessors
		bool is_open() const override;
		bool is_eof() const override;
		bool is_good() const;
		
		// file methods
		bool open() override;
		bool close() override;
			
		// reading methods
		// FileException occurs if reading failed
		// read_file method does not require to open/close the file
		bool read_line(EnhancedString& line);
		void read_file(std::vector<EnhancedString>& fileLines);
		
		// single key value file reading method
		// separator is used to separate key from the value
		// data must be of the form:
		// key[separator]value
		// separator can be repeated (eg, ::)
		void read_single_key_value_file(std::map<std::string,EnhancedString>&
                                        keys,
                                        const char& separator,
                                        const char& commentchar);

        // static methods
        static
        std::vector<std::string> read_file(const std::string& fname,
                                           const std::string& dir);
        static
        KeyValueData read_key_value_file(const std::string& fname,
                                         const std::string& dir,
                                         const char& sep,
                                         const char& comment);
		
		
	private:
		
		std::ifstream	m_input;
		
		
	};


    // key finder method for map<string,string>
    template <class KeyType> inline
    bool find_key_value(const KeyValueData& data, const std::string& key,
                        KeyType& key_value) {

        // look for key
        KeyValueData::const_iterator key_it;
        key_it = data.find(key);

        // check if the key exists and return the value otherwise return false
        if (key_it != data.end()) {
            key_value =
            EnhancedString::convert_from_string<KeyType>(key_it->second);
            return true;
        }
        else {
            key_value = KeyType();
            return false;
        };

    };
	
	
}


#endif	// DNASim_InputFileHandler_h
