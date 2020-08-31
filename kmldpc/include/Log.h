#ifndef KMLDPC_LOG_H
#define KMLDPC_LOG_H
#include <iostream>
#include <string>
#include <fstream>
#include <memory>
#include <cstring>
#include <chrono>
#include <ctime>

#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif

#define LOG(level) \
if (level > kmldpc::Log::get().getLevel()) ; \
else kmldpc::Log::get().getStream() << '[' << kmldpc::Log::get().getCurrentSystemTime() << ']' \
<< '[' << __FILENAME__ << ':' << std::dec << __LINE__ << "] "

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
        std::string getCurrentSystemTime();

        static Log& get();
    private:
        Level _logLevel;
        std::ostream* _logStream;
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
        std::streambuf* _sbtofile;
        std::streambuf* _sbtoscreen;
    };

    class TeeStream : public std::ostream
    {
    public:
        TeeStream(std::ostream& o1, std::ostream& o2);
    private:
        TeeBuf _tbuf;
    };
}


#endif //KMLDPC_LOG_H
