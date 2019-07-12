// test_EnhancedString class
// Nicolas Clauvelin


#include <EnhancedString.h>
#include <test_EnhancedString.h>
using namespace DNASim;


namespace {
	
	
	// MethodEraseBlankCharacters
	TEST_F(EnhancedStringTest, MethodEraseBlankCharacters) {
    
		const std::string s = "the lazy brown fox jumps";
		EnhancedString es(s);
		es.erase_blank_characters();
		const std::string clean_s = std::string(es);
		
		EXPECT_STREQ("thelazybrownfoxjumps", clean_s.c_str());
		
	};
	
	
	// MethodEraseCharacters
	TEST_F(EnhancedStringTest, MethodEraseCharacters) {
		
		const std::string s = "ladi,lada,,ladada";
		EnhancedString es(s);
		es.erase_character(',');
		const std::string clean_s = std::string(es);

		EXPECT_STREQ("ladiladaladada", clean_s.c_str());
		
	};
	
	
	// MethodReplaceCharacters
	TEST_F(EnhancedStringTest, MethodReplaceCharacters) {
		
		const std::string s = "ladi,lada,,ladada";
		EnhancedString es(s);
		es.replace_character(',', '-');
		const std::string clean_s = std::string(es);
		
		EXPECT_STREQ("ladi-lada--ladada", clean_s.c_str());
		
	};
	
	
	// StringTokenizingFunction
	TEST_F(EnhancedStringTest, StringTokenizingFunction) {
		
		const std::string s = "apple::pear::orange::grape";
		
		std::vector<std::string> tokens;
		bool flag = EnhancedString::tokenize_string(s, tokens, ':');
		
		EXPECT_TRUE(flag);
		EXPECT_EQ(4, (Integer)tokens.size());
		EXPECT_STREQ("apple", tokens[0].c_str());
		EXPECT_STREQ("pear", tokens[1].c_str());
		EXPECT_STREQ("orange", tokens[2].c_str());
		EXPECT_STREQ("grape", tokens[3].c_str());
		
	};
	
	
}

