#include "log.h"

namespace kmldpc
{
    Log::~Log() {}

    Log& Log::get() {
        static Log instance;
        return instance;
    }

    std::ostream& Log::getStream() {
        return *_logStream;
    }

    std::ostream& Log::getStream(bool flag) {
        auto ptr = (TeeStream* )_logStream;
        ptr->getTeeBuf().setFlag(flag);
        return *_logStream;
    }

    std::string Log::getCurrentSystemTime() {
        std::time_t secSinceEpoch = std::chrono::system_clock::to_time_t(
                std::chrono::system_clock::now());
        struct tm* calendarTime = localtime(&secSinceEpoch);
        char usrdefFormat[50] = { 0 };
        strftime(usrdefFormat, 50,
                "%Y-%m-%d %H:%M:%S", calendarTime);
        return std::string(usrdefFormat);
    }

    void Log::setLogStream(std::ostream &stream) {
        this->_logStream = &stream;
    }

    Log& Log::setLevel(Level level) {
        _logLevel = level;
        return *this;
    }

    Level Log::getLevel() {
        return _logLevel;
    }

    TeeBuf::TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen) :
            _sbtofile(sbtofile),
            _sbtoscreen(sbtoscreen),
            _flag(true)
    {}

    void TeeBuf::setFlag(bool flag) {
        this->_flag = flag;
    }

    int TeeBuf::overflow(int c) {
        if (c == EOF) {
            return !EOF;
        } else {
            const int r1 = _sbtofile->sputc(c);
            if (_flag) {
                const int r2 = _sbtoscreen->sputc(c);
                return r1 == EOF || r2 == EOF ? EOF : c;
            } else {
                return r1 == EOF ? EOF : c;
            }
        }
    }

    int TeeBuf::sync() {
        const int r1 = _sbtofile->pubsync();
        if (_flag) {
            const int r2 = _sbtoscreen->pubsync();
            return r1 == 0 && r2 == 0 ? 0 : -1;
        } else {
            return r1 == 0 ? 0 : -1;
        }
    }

    TeeStream::TeeStream(std::ostream &o1, std::ostream &o2) :
            std::ostream(&_tbuf),
            _tbuf(o1.rdbuf(), o2.rdbuf())
    {}

    TeeBuf& TeeStream::getTeeBuf() {
        return _tbuf;
    }
}