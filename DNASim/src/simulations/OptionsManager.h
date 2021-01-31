// OptionsManager class
// Nicolas Clauvelin


// simple option manager implementation

// the options are stored as pairs of strings and convert to the desired type
// when requested

// for boolean options, set the string to "1" for true and "0" for false


#ifndef DNASim_OptionsManager_h
#define DNASim_OptionsManager_h


#include "file_io/EnhancedString.h"
#include "file_io/InputFileHandler.h"


namespace DNASim {
	
	
	class OptionsManager {
		
		
	public:
		
		// constructors
		OptionsManager() = default;
        OptionsManager(const std::map<std::string,std::string>& options) :
        m_options(options) {};
		~OptionsManager() = default;
		
		// number of options accessor
		inline Size n_options() const {
			return m_options.size();
		};
		
		// set options list method
		inline void set_options_list(const std::vector<std::string>& options) {
			m_options.clear();
			for (auto it = options.begin(), end = options.end();
                 it != end;
                 ++it)
				m_options[*it] = std::string();
		};
		
		// options setting method
		inline bool set_option(const std::string& option_name,
							   const std::string& option_value) {
			
			// look for option
			auto it = m_options.find(option_name);
			
			// set option if it exists otherwise return false
			if (it != m_options.end()) {
                it->second = option_value;
                return true;
            };

            return false;

		};
		
		// options getting method
		template <typename ValueType> inline
		bool get_option(const std::string& option_name, ValueType& value) const
        {

            // lookup
            bool lkup_res = find_key_value(m_options, option_name, value);

            // return false if the option is not found
			return lkup_res;
			
		};
		
		// iterator over options
		std::map<std::string,std::string>::const_iterator begin() const {
			return m_options.begin();
		};
		std::map<std::string,std::string>::const_iterator end() const {
			return m_options.end();
		};
		
		
	private:
		
		KeyValueData m_options;
		
		
	};
			

}


#endif	// DNASim_OptionsManager_h
