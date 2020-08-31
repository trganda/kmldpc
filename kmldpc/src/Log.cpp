#include "Log.h"

namespace kmldpc
{
    Log::~Log() {}

    Log& Log::get() {
        static Log instance;
        return instance;
    }

    std::ostream& Log::getStream() {
        return *m_logStream;
    }

    void Log::setLogStream(std::ostream &stream) {
        this->m_logStream = &stream;
    }

    Log& Log::setLevel(Level level) {
        m_logLevel = level;
        return *this;
    }

    Level Log::getLevel() {
        return m_logLevel;
    }

    TeeBuf::TeeBuf(std::streambuf *sbtofile, std::streambuf *sbtoscreen) :
        m_sbtofile(sbtofile),
        m_sbtoscreen(sbtoscreen)
    {}

    int TeeBuf::overflow(int c) {
        if (c == EOF) {
            return !EOF;
        } else {
            const int r1 = m_sbtofile->sputc(c);
            const int r2 = m_sbtoscreen->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    int TeeBuf::sync() {
        const int r1 = m_sbtofile->pubsync();
        const int r2 = m_sbtoscreen->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }

    TeeStream::TeeStream(std::ostream &o1, std::ostream &o2) :
        std::ostream(&m_tbuf),
        m_tbuf(o1.rdbuf(), o2.rdbuf())
    {}

}