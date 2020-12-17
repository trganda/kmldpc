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
            explicit TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen)
                    : strbuf_to_file_(sbtofile),
                      strbuf_to_stdout_(sbtoscreen),
                      _flag(true) {}

            void SetFlag(bool flag) {
                this->_flag = flag;
            }

        private:
            // This tee buffer has no buffer. So every character "overflows"
            // and can be put directly into the teed buffers.
            int overflow(int c) override {
                if (c == EOF) {
                    return !EOF;
                } else {
                    const int r1 = strbuf_to_file_->sputc(c);
                    if (_flag) {
                        const int r2 = strbuf_to_stdout_->sputc(c);
                        return r1 == EOF || r2 == EOF ? EOF : c;
                    } else {
                        return r1 == EOF ? EOF : c;
                    }
                }
            }

            // Sync both teed buffers.
            int sync() override {
                const int r1 = strbuf_to_file_->pubsync();
                if (_flag) {
                    const int r2 = strbuf_to_stdout_->pubsync();
                    return r1 == 0 && r2 == 0 ? 0 : -1;
                } else {
                    return r1 == 0 ? 0 : -1;
                }
            }

        private:
            std::streambuf *strbuf_to_file_;
            std::streambuf *strbuf_to_stdout_;
            bool _flag;
        };

        class TeeStream : public std::ostream {
        public:
            explicit TeeStream(std::ostream &o1, std::ostream &o2)
                    : std::ostream(&tbuf_),
                      tbuf_(o1.rdbuf(), o2.rdbuf()) {}

            TeeBuf &GetTeeBuf() {
                return tbuf_;
            }

        private:
            TeeBuf tbuf_;
        };

        class Log {
        public:
            ~Log() = default;

            void SetLogStream(std::ostream &stream) {
                this->log_stream_ = &stream;
            }

            Log &SetLevel(Level level) {
                log_level_ = level;
                return *this;
            }

            Level GetLevel() {
                return log_level_;
            }

            std::ostream &GetStream() {
                return *log_stream_;
            }

            std::ostream &GetStream(bool flag) {
                auto ptr = (TeeStream *) log_stream_;
                ptr->GetTeeBuf().SetFlag(flag);
                return *log_stream_;
            }

            static std::string GetCurrentSystemTime() {
                std::time_t secSinceEpoch = std::chrono::system_clock::to_time_t(
                        std::chrono::system_clock::now());
                struct tm *calendarTime = localtime(&secSinceEpoch);
                char usrdefFormat[50] = {0};
                strftime(usrdefFormat, 50,
                         "%Y-%m-%d %H:%M:%S", calendarTime);
                return std::string(usrdefFormat);
            }

            static Log &Get() {
                static Log instance;
                return instance;
            }

        private:
            Level log_level_;
            std::ostream *log_stream_;
        };

    }

}   // namespace


#endif //KMLDPC_LOG_H
