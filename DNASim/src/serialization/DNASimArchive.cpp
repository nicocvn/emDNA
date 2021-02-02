// DNASimArchive class
// Nicolas Clauvelin


#include "serialization/DNASimArchive.h"


namespace DNASim {


    // output archive initialization method - binary
    template <>
    void OutputArchive<CerealBinaryOutput>::init_archive(const bool&
                                                         append_flag) {

        // output stream
        if (!append_flag)
            m_archive_stream_ptr =
            std::unique_ptr<std::ofstream>(new
                                           std::ofstream(m_archive_filename,
                                                         std::
                                                         ios_base::binary));
        else
            m_archive_stream_ptr =
            std::unique_ptr<std::ofstream>(new
                                           std::ofstream(m_archive_filename,
                                                         std::ios_base::binary |
                                                         std::ios::app));

        // archive
        m_archive_ptr =
        std::
        unique_ptr<CerealBinaryOutput>(new
                                       CerealBinaryOutput(*m_archive_stream_ptr)
                                       );

    };


    // output archive initialization method - JSON
    template <>
    void OutputArchive<CerealJSONOutput>::init_archive(const bool& append_flag)
    {

        // output stream
        if (!append_flag)
            m_archive_stream_ptr =
            std::unique_ptr<std::ofstream>(new
                                           std::ofstream(m_archive_filename));
        else
            m_archive_stream_ptr =
            std::unique_ptr<std::ofstream>(new std::ofstream(m_archive_filename,
                                                             std::ios::app));

        // archive
        m_archive_ptr =
        std::unique_ptr<CerealJSONOutput>
        (new CerealJSONOutput(*m_archive_stream_ptr));

    };


    // output archive initialization method - XML
    template <>
    void OutputArchive<CerealXMLOutput>::init_archive(const bool& append_flag)
    {

        // output stream
        if (!append_flag)
            m_archive_stream_ptr =
            std::unique_ptr<std::ofstream>(new
                                           std::ofstream(m_archive_filename));
        else
            m_archive_stream_ptr =
            std::unique_ptr<std::ofstream>(new std::ofstream(m_archive_filename,
                                                             std::ios::app));

        // archive
        m_archive_ptr =
        std::unique_ptr<CerealXMLOutput>(new
                                         CerealXMLOutput(*m_archive_stream_ptr)
                                         );
        
    };


    // input archive initialization method - binary
    template <>
    void InputArchive<CerealBinaryInput>::init_archive() {

        // input stream
        m_archive_stream_ptr =
        std::unique_ptr<std::ifstream>(new
                                       std::ifstream(m_archive_filename,
                                                     std::ios_base::binary));

        // archive
        m_archive_ptr =
        std::unique_ptr<CerealBinaryInput>
        (new CerealBinaryInput(*m_archive_stream_ptr));

    };


    // input archive initialization method - JSON
    template <>
    void InputArchive<CerealJSONInput>::init_archive() {

        // input stream
        m_archive_stream_ptr =
        std::unique_ptr<std::ifstream>(new std::ifstream(m_archive_filename));

        // archive
        m_archive_ptr =
        std::unique_ptr<CerealJSONInput>(new
                                         CerealJSONInput(*m_archive_stream_ptr)
                                         );
    };


    // input archive initialization method - XML
    template <>
    void InputArchive<CerealXMLInput>::init_archive() {

        // input stream
        m_archive_stream_ptr =
        std::unique_ptr<std::ifstream>(new std::ifstream(m_archive_filename));

        // archive
        m_archive_ptr =
        std::unique_ptr<CerealXMLInput>(new
                                        CerealXMLInput(*m_archive_stream_ptr));
    };


}
