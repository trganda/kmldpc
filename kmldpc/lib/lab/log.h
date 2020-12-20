#ifndef LAB_LOG_HPP
#define LAB_LOG_HPP

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

#define LOG(level, flag) \
if (level > lab::logger::Log::Get().GetLevel()) ; \
else lab::logger::Log::Get().GetStream(flag) << '[' << lab::logger::Log::GetCurrentSystemTime() << ']' \
<< '[' << __FILENAME__ << ':' << std::dec << __LINE__ << "] "

namespace lab {

    namespace logger {

        enum Level {
            None,
            Error,
            Info,
            InforVerbose,
        };

        //Courtesy of http://wordaligned.org/articles/cpp-streambufs#toctee-streams
        class TeeBuf : public std::streambuf {
        public:
            // Construct a streambuf which tees output to both input
            // streambufs.
            explicit TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen);
            void SetFlag(bool flag);

        private:
            // This tee buffer has no buffer. So every character "overflows"
            // and can be put directly into the teed buffers.
            int overflow(int c) override;

            // Sync both teed buffers.
            int sync() override;

        private:
            std::streambuf *strbuf_to_file_;
            std::streambuf *strbuf_to_stdout_;
            bool _flag;
        };

        class TeeStream : public std::ostream {
        public:
            explicit TeeStream(std::ostream &o1, std::ostream &o2);
            TeeBuf &GetTeeBuf();

        private:
            TeeBuf tbuf_;
        };

        class Log {
        public:
            ~Log() = default;

            void SetLogStream(std::ostream &stream);
            Log &SetLevel(Level level);
            Level GetLevel();
            std::ostream &GetStream();
            std::ostream &GetStream(bool flag);
            static std::string GetCurrentSystemTime();
            static Log &Get();

        private:
            Level log_level_;
            std::ostream *log_stream_;
        };

    }

}   // namespace


#endif //KMLDPC_LOG_H
