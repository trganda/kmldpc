#ifndef KMLDPC_LOG_H
#define KMLDPC_LOG_H
#include <iostream>
#include <string>
#include <fstream>
#include <memory>
#include <cstring>

#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif

#define LOG(level) \
if (level > kmldpc::Log::get().getLevel()) ; \
else kmldpc::Log::get().getStream() << '[' << __FILENAME__ << ":" << std::dec << __LINE__ << "] "

namespace kmldpc
{
    enum Level
    {
        None,
        Error,
        Info,
        InforVerbose,
    };

    class Log
    {
    public:
        ~Log();
        void setLogStream(std::ostream& stream);
        Log& setLevel(Level level);
        Level getLevel();

        std::ostream& getStream();

        static Log& get();
    private:
        Level m_logLevel;
        std::ostream* m_logStream;
    };
    //Courtesy of http://wordaligned.org/articles/cpp-streambufs#toctee-streams
    class TeeBuf : public std::streambuf
    {
    public:
        // Construct a streambuf which tees output to both input
        // streambufs.
        TeeBuf(std::streambuf* sbtofile, std::streambuf* sbtoscreen);
    private:
        // This tee buffer has no buffer. So every character "overflows"
        // and can be put directly into the teed buffers.
        virtual int overflow(int c);
        // Sync both teed buffers.
        virtual int sync();
    private:
        std::streambuf* m_sbtofile;
        std::streambuf* m_sbtoscreen;
    };

    class TeeStream : public std::ostream
    {
    public:
        TeeStream(std::ostream& o1, std::ostream& o2);
    private:
        TeeBuf m_tbuf;
    };
}


#endif //KMLDPC_LOG_H
