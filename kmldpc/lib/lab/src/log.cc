#include "log.h"

namespace lab {
namespace logger {
TeeBuf::TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen)
    : strbuf_to_file_(sbtofile),
      strbuf_to_stdout_(sbtoscreen),
      _flag(true) {}

void TeeBuf::SetFlag(bool flag) {
  this->_flag = flag;
}

int TeeBuf::overflow(int c) {
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

int TeeBuf::sync() {
  const int r1 = strbuf_to_file_->pubsync();
  if (_flag) {
    const int r2 = strbuf_to_stdout_->pubsync();
    return r1 == 0 && r2 == 0 ? 0 : -1;
  } else {
    return r1 == 0 ? 0 : -1;
  }
}

TeeStream::TeeStream(std::ostream &o1, std::ostream &o2)
    : std::ostream(&tbuf_),
      tbuf_(o1.rdbuf(), o2.rdbuf()) {}

TeeBuf &TeeStream::GetTeeBuf() {
  return tbuf_;
}

void Log::SetLogStream(std::ostream &stream) {
  this->log_stream_ = &stream;
}

Log &Log::SetLevel(Level level) {
  log_level_ = level;
  return *this;
}

Level Log::GetLevel() {
  return log_level_;
}

std::ostream &Log::GetStream() {
  return *log_stream_;
}

std::ostream &Log::GetStream(bool flag) {
  auto ptr = (TeeStream *) log_stream_;
  ptr->GetTeeBuf().SetFlag(flag);
  return *log_stream_;
}

std::string Log::GetCurrentSystemTime() {
  std::time_t secSinceEpoch = std::chrono::system_clock::to_time_t(
      std::chrono::system_clock::now());
  struct tm *calendarTime = localtime(&secSinceEpoch);
  char usrdefFormat[50] = {0};
  strftime(usrdefFormat, 50,
           "%Y-%m-%d %H:%M:%S", calendarTime);
  return std::string(usrdefFormat);
}

Log &Log::Get() {
  static Log instance;
  return instance;
}
}
}   // namespace
