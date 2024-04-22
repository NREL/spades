#include "Utilities.H"

namespace spades {
// skip to next line in Header
void goto_next_line(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

// equivalent to amrex::ParallelDescriptor::ReadAndBcastFile but for all ranks
void ReadFile(
    const std::string& filename,
    amrex::Vector<char>& charBuf,
    bool bExitOnError)
{
    enum { IO_Buffer_Size = 262144 * 8 };

#ifdef BL_SETBUF_SIGNED_CHAR
    using Setbuf_Char_Type = signed char;
#else
    using Setbuf_Char_Type = char;
#endif

    amrex::Vector<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    amrex::Long fileLength(0), fileLengthPadded(0);

    std::ifstream iss;

    iss.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    iss.open(filename.c_str(), std::ios::in);
    if (!iss.good()) {
        if (bExitOnError) {
            amrex::FileOpenFailed(filename);
        } else {
            fileLength = -1;
        }
    } else {
        iss.seekg(0, std::ios::end);
        fileLength = static_cast<std::streamoff>(iss.tellg());
        iss.seekg(0, std::ios::beg);
    }

    if (fileLength == -1) {
        return;
    }

    fileLengthPadded = fileLength + 1;
    //    fileLengthPadded += fileLengthPadded % 8;
    charBuf.resize(fileLengthPadded);

    iss.read(charBuf.dataPtr(), fileLength);
    iss.close();

    charBuf[fileLength] = '\0';
}
} // namespace spades
