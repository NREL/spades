#include "Utilities.H"

namespace spades {
// skip to next line in Header
void goto_next_line(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

// equivalent to amrex::ParallelDescriptor::ReadAndBcastFile but for all ranks
void read_file(
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

    amrex::Long file_length(0), file_length_padded(0);

    std::ifstream iss;

    iss.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    iss.open(filename.c_str(), std::ios::in);
    if (!iss.good()) {
        if (bExitOnError) {
            amrex::FileOpenFailed(filename);
        } else {
            file_length = -1;
        }
    } else {
        iss.seekg(0, std::ios::end);
        file_length = static_cast<std::streamoff>(iss.tellg());
        iss.seekg(0, std::ios::beg);
    }

    if (file_length == -1) {
        return;
    }

    file_length_padded = file_length + 1;
    //    file_length_padded += file_length_padded % 8;
    charBuf.resize(file_length_padded);

    iss.read(charBuf.dataPtr(), file_length);
    iss.close();

    charBuf[file_length] = '\0';
}
} // namespace spades
