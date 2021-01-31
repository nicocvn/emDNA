// EnhancedString class
// Nicolas Clauvelin


#include "file_io/EnhancedString.h"


namespace DNASim {
	
	
	// empty string checking method
	bool EnhancedString::is_empty() const {
		if (m_s.empty() == 0)
			return true;
		return false;
	};
	
	
	// size of the string accessor
	UInteger EnhancedString::size() const {
		return (UInteger)m_s.size();
	};
	
	
	// std::string cast operator
	EnhancedString::operator std::string() const {
		return m_s;
	};
	
	
	// UInteger cast operator
	EnhancedString::operator UInteger() const {
		return convert<UInteger>();
	};
	
	
	// Integer cast operator
	EnhancedString::operator Integer() const {
		return convert<Integer>();
	};
	
	
	// Real cast operator
	EnhancedString::operator Real() const {
		return convert<Real>();
	};
	
	
	// bool cast operator
	EnhancedString::operator bool() const {
		if (m_s=="true")
			return true;
		else
			return false;
	};
	
	
	// std::string operator =
	// copy the content of s in the enhanced string (exact copy)
	EnhancedString&	EnhancedString::operator=(const std::string& s) {
		m_s = s;
		return *this;
	};
	
	
	// std::string operator +=
	// concatenate s with the enhanced string
	EnhancedString&	EnhancedString::operator+=(const std::string& s) {
		m_s += s;
		return *this;
	};
	
	
	// std::string operator <<
	// same as operator =
	EnhancedString&	EnhancedString::operator<<(const std::string& s) {
		m_s = s;
		return *this;
	};

	
	// operator +=
	// string concatenation
	EnhancedString&	EnhancedString::operator+=(const EnhancedString& s) {
		m_s += s.m_s;
		return *this;
	};
	
	
	// << operator
	// same as operator =
	EnhancedString&	EnhancedString::operator<<(const EnhancedString& s) {
		m_s = s.m_s;
		return *this;
	};
	
	
	// == operator
	// equality operator
	bool EnhancedString::operator==(const EnhancedString& s) {
		if (m_s == s.m_s)
			return true;
		else
			return false;
	};
	
	
	// tokenizing method
	// if no tokens are found the full string is put in std::vector tokens and the
	// method returns false
	bool EnhancedString::tokenize(std::vector<EnhancedString>& tokens,
								  const char& delimiter) const {
		
		// clean tokens
		tokens.clear();
		
		// skip delimiters at begining
		std::string::size_type lastpos = m_s.find_first_not_of(delimiter, 0);
		
		// find first non-delimiter digit
		std::string::size_type pos = m_s.find_first_of(delimiter, lastpos);
		
		// tokenizing loop
		while (std::string::npos != pos || std::string::npos != lastpos) {
			
			// add token to the vector
			tokens.push_back(EnhancedString(m_s.substr(lastpos,
													   pos - lastpos)));
			
			// skip delimiters
			lastpos = m_s.find_first_not_of(delimiter, pos);
			
			// find next non-delimiter digit
			pos = m_s.find_first_of(delimiter, lastpos);
			
		};
		
		// no tokenizing
		if (tokens.size() == 1)
			return false;
		
		return true;
		
	};
	
	
	// std::string tokenizing method
	// if no tokens are found the full string is put in std::vector tokens and
	// the method returns false
	bool EnhancedString::tokenize(std::vector<std::string>& tokens,
								  const char& delimiter) const {
		
		// clean tokens
		tokens.clear();
		
		// skip delimiters at begining
		std::string::size_type lastpos = m_s.find_first_not_of(delimiter, 0);
		
		// find first non-delimiter digit
		std::string::size_type pos = m_s.find_first_of(delimiter, lastpos);
		
		// tokenizing loop
		while (std::string::npos != pos || std::string::npos != lastpos) {
			
			// add token to the vector
			tokens.push_back(m_s.substr(lastpos, pos - lastpos));
			
			// skip delimiters
			lastpos = m_s.find_first_not_of(delimiter, pos);
			
			// find next non-delimiter digit
			pos = m_s.find_first_of(delimiter, lastpos);
			
		};
		
		// no tokenizing
		if (tokens.size() == 1)
			return false;
		
		return true;
		
	};
	
	
	// blank characters erase method
	// delete white spaces and tabs from the string
	void EnhancedString::erase_blank_characters() {
		
		// use a condensed form with remove_if and predicates
		m_s.erase(std::remove_if(m_s.begin(), m_s.end(), isblank), m_s.end());
		
	};
	
	
	// erase character method
	// erase all the occurrences of c in the enhanced string
	void EnhancedString::erase_character(const char& c) {

        // find first occurence
        Size found = m_s.find(c);
        
        // loop over the string and delete all occurences
        while (found != std::string::npos) {
            m_s.erase(found, 1);
            found = m_s.find(c);
        };
		
	};
    void EnhancedString::
    erase_characters(const std::vector<char>& cvec) {
        for (auto c : cvec)
            erase_character(c);
    };
	
	
	// replace character method
	// replace all the occurrences of c in the enhanced string
	void EnhancedString::replace_character(const char& c,
                                           const char& newc) {
		
		std::replace(m_s.begin(), m_s.end(), c, newc);
		
	};


	// friend concatenation operator
	EnhancedString operator+(const EnhancedString& s1,
							 const EnhancedString& s2) {

		EnhancedString es(s1);
		es += s2;

		return es;

	};


	// output operator
	std::ostream& operator<<(std::ostream& os, const EnhancedString& s) {

		os << s.m_s;

		return os;

	};
	
	
	// static std::string tokenizing method
	bool EnhancedString::tokenize_string(const std::string& s,
										 std::vector<std::string>& tokens,
										 const char& delimiter) {
		
		EnhancedString es(s);
		return es.tokenize(tokens, delimiter);
		
	};


    // static std::string erase character method
    std::string EnhancedString::erase_character(const std::string& s,
                                                const char& c) {
        EnhancedString es(s);
        es.erase_character(c);
        return std::string(es);
    };
    
    
    // friend whitespace trimming method
    std::string EnhancedString::remove_spaces(const std::string& s) {
        EnhancedString es(s);
        es.erase_blank_characters();
        return es.m_s;
    };
    
    
    // friend character erasing method
    std::string EnhancedString::remove_character(const std::string& s,
                                                 const char& c) {
        EnhancedString es(s);
        es.erase_character(c);
        return es.m_s;
    };


    // friend characters removing method
    std::string EnhancedString::remove_characters(const std::string& s,
                                                  const std::vector<char>&
                                                  chars) {
        EnhancedString es(s);
        es.erase_characters(chars);
        return es.m_s;
    };


    // friend string tokenizing method
    std::vector<std::string>
    EnhancedString::tokenize_string(const std::string& s,
                                    const char& sep) {
        std::vector<std::string> tokens;
        EnhancedString::tokenize_string(s, tokens, sep);
        return tokens;
    };

		
}
