// DNASim Debug header
// Nicolas Clauvelin


// project-wide header defining macros used for assertion and debugging


#ifndef DNASim_DNASim_Debug_h
#define DNASim_DNASim_Debug_h


// exit exception structure
struct DNASim_ExitException {
    int _exit_code;
    DNASim_ExitException(int code) : _exit_code(code) {};
};


namespace {

	// DNASim assert macro
	#define DS_ASSERT(_expr, _message) \
	do { \
		if (_expr != true) {\
        std::cerr << "\n"; \
        std::cerr << "[[ ASSERT ERROR ]]" << "\n"; \
		std::cerr << "assert: " << #_expr << "\n"; \
		std::cerr << _message; \
		std::cerr << "\n"; \
		std::cerr << "from: " <<__PRETTY_FUNCTION__ \
		<< " @line " <<__LINE__ << "\n"; \
		throw DNASim_ExitException(-1); \
	} else { break; }} while(0);
	
	
}


#endif	// DNASim_DNASim_Debug_h
