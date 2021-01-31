// DNASimArchive class
// Nicolas Clauvelin


// cereal archive wrapper class for binary, JSON, XML archives

// the archive stream is implemented as a unique_ptr to avoid issues with
// GCC loosy implementation of iostream move operators


#ifndef DNASim_DNASimArchive_h
#define DNASim_DNASimArchive_h


#include "DNASim_Includes.h"
#include "serialization/Cereal_Includes.h"


namespace DNASim {


    // template output archive class
    template <class ArchiveType>
    class OutputArchive {


    public:

        // constructor
        OutputArchive(const std::string& filename,
                      const bool& append);
        ~OutputArchive();

        // archive file name accessor
        const std::string& archive_filename() const;

        // save methods
        template <class ObjectType> void save(const ObjectType& object) {
            (*m_archive_ptr)(object);
        };
        template <class PtrType>
        void save(const std::shared_ptr<PtrType>& object_ptr) {
            (*m_archive_ptr)(object_ptr);
        };


    private:

        // deleted constructors and operators
        OutputArchive() = delete;
        OutputArchive(const OutputArchive<ArchiveType>& arxiv) = delete;
        OutputArchive(OutputArchive<ArchiveType>&& arxiv) = delete;
        OutputArchive&
        operator=(const OutputArchive<ArchiveType>& arxiv) = delete;
        OutputArchive& operator=(OutputArchive<ArchiveType>&&
                                 arxiv) = delete;

        // initialize archive method
        void init_archive(const bool& append_flag);

        // archive file name and stream
        std::string m_archive_filename;
        std::unique_ptr<std::ofstream> m_archive_stream_ptr;

        // cereal binary archive
        std::unique_ptr<ArchiveType> m_archive_ptr;


    };


    // template input archive
    template <class ArchiveType>
    class InputArchive {


    public:

        // constructor
        InputArchive(const std::string& filename);
        ~InputArchive();

        // archive file name accessor
        const std::string& archive_filename() const;

        // load method
        template <class ObjectType> void load(ObjectType& object) {
            (*m_archive_ptr)(object);
        };
        template <class PtrType>
        void load(std::shared_ptr<PtrType>& object_ptr) {
            (*m_archive_ptr)(object_ptr);
        };


    private:

        // deleted constructors and operators
        InputArchive() = delete;
        InputArchive(const InputArchive<ArchiveType>& arxiv) = delete;
        InputArchive(InputArchive<ArchiveType>&& arxiv) = delete;
        InputArchive&
        operator=(const InputArchive<ArchiveType>& arxiv) = delete;
        InputArchive& operator=(InputArchive<ArchiveType>&& arxiv) = delete;

        // initialize archive method
        void init_archive();

        // archive file name and stream
        std::string m_archive_filename;
        std::unique_ptr<std::ifstream> m_archive_stream_ptr;

        // cereal binary archive
        std::unique_ptr<ArchiveType> m_archive_ptr;


    };


    // class constructor with archive initialization
    template <class ArchiveType>
    OutputArchive<ArchiveType>::OutputArchive(const std::string& filename,
                                              const bool& append) :
    m_archive_filename(filename),
    m_archive_stream_ptr(nullptr),
    m_archive_ptr(nullptr) {
        init_archive(append);
    };


    // class destructor
    template <class ArchiveType>
    OutputArchive<ArchiveType>::~OutputArchive() {};


    // archive file name accessor/modifier
    template <class ArchiveType>
    const std::string& OutputArchive<ArchiveType>::archive_filename() const {
        return m_archive_filename;
    };


    // class constructor with archive initialization
    template <class ArchiveType>
    InputArchive<ArchiveType>::InputArchive(const std::string& filename) :
    m_archive_filename(filename),
    m_archive_stream_ptr(nullptr),
    m_archive_ptr(nullptr) {
        init_archive();
    };


    // class destructor
    template <class ArchiveType>
    InputArchive<ArchiveType>::~InputArchive() {};


    // archive file name accessor/modifier
    template <class ArchiveType>
    const std::string& InputArchive<ArchiveType>::archive_filename() const {
        return m_archive_filename;
    };


}


#endif  // DNASim_DNASimArchive_h
