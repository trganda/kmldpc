#ifndef LAB_LOG_H
#define LAB_LOG_H

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
if (level > lab::logger::Log::get().log_level()) ; \
else lab::logger::Log::get().log_stream(flag) << '[' << lab::logger::Log::get_time() << ']' \
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
    void set_flag(bool flag);
 private:
    // This tee buffer has no buffer. So every character "overflows"
    // and can be put directly into the teed buffers.
    int overflow(int c) override;
    // Sync both teed buffers.
    int sync() override;
 private:
    std::streambuf *strbuf_to_file_;
    std::streambuf *strbuf_to_stdout_;
    bool flag_;
};
class TeeStream : public std::ostream {
 public:
    explicit TeeStream(std::ostream &o1, std::ostream &o2);
    TeeBuf &tbuf();
 private:
    TeeBuf tbuf_;
};
class Log {
 public:
    ~Log() = default;
    void set_log_stream(std::ostream &stream);
    Log &set_log_level(Level level);
    Level log_level() const;
    std::ostream &log_stream() const;
    std::ostream &log_stream(bool flag) const;
    static Log &get();
    static std::string get_time();
 private:
    Level log_level_;
    std::ostream *log_stream_;
};
}
} // namespace lab
#endif
