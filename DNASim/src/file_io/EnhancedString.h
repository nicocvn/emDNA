// EnhancedString class
// Nicolas Clauvelin


// enhanced string class
// provides tokenizing, conversion and erase/replace characters methods


#ifndef DNASim_Enhanced_String_h
#define DNASim_Enhanced_String_h


#include "DNASim_Includes.h"


namespace DNASim {
	
	
	class EnhancedString {
		
		
	public:
		
		// constructors
		// default constructor creates an empty string
		EnhancedString() = default;
		explicit EnhancedString(const std::string& s) : m_s(s) {};
        explicit EnhancedString(std::string&& s) : m_s(std::move(s)) {};
        EnhancedString(const EnhancedString& es) = default;
        EnhancedString(EnhancedString&& es) = default;
		~EnhancedString() = default;
		
		// string properties
		bool is_empty() const;
		UInteger size() const;
		
		// cast operators
		explicit operator std::string() const;
		explicit operator UInteger() const;
		explicit operator Integer() const;
		explicit operator Real() const;
		explicit operator bool() const;
		
		// operators for std::string
		EnhancedString&	operator=(const std::string& s);
		EnhancedString&	operator+=(const std::string& s);
		EnhancedString&	operator<<(const std::string& s);
		
		// operators
		EnhancedString&	operator=(const EnhancedString& s) = default;
        EnhancedString&	operator=(EnhancedString&& s) = default;
		EnhancedString&	operator+=(const EnhancedString& s);
		EnhancedString&	operator<<(const EnhancedString& s);
		bool operator==(const EnhancedString& s);
		
		// string conversion method
		template <class T> T convert() const {
			std::istringstream stream(m_s);
			T t;
			stream >> t;
			return t;
		};
		
		// string tokenizing methods
		// tokenizing is implemented to return tokens as either EnhancedString
		// or regular std::string
		// return false if no tokens are found
		bool tokenize(std::vector<EnhancedString>& tokens,
                      const char& delimiter) const;
		bool tokenize(std::vector<std::string>& tokens,
                      const char& delimiter) const;
		
		// string manipulation methods
		void erase_blank_characters();
		void erase_character(const char& c);
        void erase_characters(const std::vector<char>& cvec);
		void replace_character(const char& c, const char& newc);

        // friend operator for string concatenation
		friend EnhancedString operator+(const EnhancedString& s1,
										const EnhancedString& s2);

		// output operator
		friend std::ostream& operator<<(std::ostream& os,
										const EnhancedString& s);
		
		// static methods for std::string
		static bool tokenize_string(const std::string& s,
									std::vector<std::string>& tokens,
									const char& delimiter);
        static std::string erase_character(const std::string& s,
                                           const char& c);

        // static methods for string conversion
        template <class T>
        static T convert_from_string(const std::string& T_str);
        template <class T>
        static std::string convert_to_string(const T& t);
        
        // static methods for string manipulation
        static std::string remove_spaces(const std::string& s);
        static std::string remove_character(const std::string& s,
                                            const char& c);
        static std::string remove_characters(const std::string& s,
                                             const std::vector<char>& chars);
        static std::vector<std::string> tokenize_string(const std::string& s,
                                                        const char& sep);
		
		
	private:
		
		std::string m_s;
		
		
	};


    // static methods for string conversion
    template <class T>
    T EnhancedString::convert_from_string(const std::string& T_str) {
        return EnhancedString(T_str).convert<T>();
    };
    template <class T>
    std::string EnhancedString::convert_to_string(const T& t) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(REAL_WIDTH) << t;
        return std::string(ss.str());
    };


    // static method for vector list parsing
    // T is required to have a std::string constructor
    template <class T>
    std::vector<T> read_vector(const std::string& data_str,
                               const char& opening,
                               const char& closing,
                               const char& separator) {

        // string copy
        std::string input_str(data_str);

        // remove first and last characters
        input_str.erase(input_str.begin());
        input_str.pop_back();

        // tokenizing
        std::vector<std::string> tokens =
        EnhancedString::tokenize_string(input_str, separator);

        // parsing
        std::vector<T> objects;
        objects.reserve(tokens.size());
        for (const std::string& tk : tokens)
            objects.push_back(T(tk));

        return objects;
        
    };

	
}


#endif	// DNASim_Enhanced_String_h
