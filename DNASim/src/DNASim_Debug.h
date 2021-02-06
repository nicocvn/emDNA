// DNASim Debug header
// Nicolas Clauvelin


// project-wide header defining macros used for assertion and debugging
#ifndef DNASim_DNASim_Debug_h
#define DNASim_DNASim_Debug_h


// exit exception structure
struct DNASim_ExitException {
    int _exit_code;
    DNASim_ExitException(int code) : _exit_code(code){};
};


namespace {

#if defined(__GNUC__) || defined(__llvm__)
#    define FUNC_MACRO __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#    define FUNC_MACRO __FUNCSIG__
#else
#    define FUNC_MACRO "????"
#endif

// DNASim assert macro
#define DS_ASSERT(_expr, _message)                                             \
    do {                                                                       \
        if (_expr != true) {                                                   \
            std::cerr << "\n";                                                 \
            std::cerr << "[[ ASSERT ERROR ]]"                                  \
                      << "\n";                                                 \
            std::cerr << "assert: " << #_expr << "\n";                         \
            std::cerr << _message;                                             \
            std::cerr << "\n";                                                 \
            std::cerr << "from: " << FUNC_MACRO << " @line " << __LINE__       \
                      << "\n";                                                 \
            throw DNASim_ExitException(-1);                                    \
        }                                                                      \
        else {                                                                 \
            break;                                                             \
        }                                                                      \
    } while (0);

}    // namespace


#endif    // DNASim_DNASim_Debug_h
