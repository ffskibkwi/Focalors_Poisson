#include "profiler.h"

const std::string cdel_str = "__cdecl ";

void ProfilerSingleton::BeginSession(const std::string& name, const std::string& filepath)
{
    if (is_writing_enabled == false)
        return;

    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_current_session)
    {
        std::cout << "Session " << std::quoted(m_current_session->Name) << " already open." << std::endl;

        InternalEndSession();
    }
    m_output_stream.open(filepath);

    if (m_output_stream.is_open())
    {
        m_current_session = new ProfilerSession({name});
        WriteHeader();
    }
    else
    {
        std::cout << "Could not open results file " << std::quoted(filepath) << "." << std::endl;
    }
}

void ProfilerSingleton::EndSession()
{
    if (is_writing_enabled == false)
        return;

    std::unique_lock<std::mutex> lock(m_mutex);
    InternalEndSession();
}

void ProfilerSingleton::Upload(const ProfileResult& result)
{
    if (is_writing_enabled == false)
        return;

    std::stringstream json;

    json << std::setprecision(3) << std::fixed;
    json << ",\n\t\t{";
    json << "\n\t\t\t\"cat\":\"function\",";
    json << "\n\t\t\t\"dur\":" << (result.ElapsedTime.count()) << ',';
    json << "\n\t\t\t\"name\":\"" << result.Name << "\",";
    json << "\n\t\t\t\"ph\":\"X\",";
    json << "\n\t\t\t\"pid\":" << mpi_rank << ",";
    json << "\n\t\t\t\"tid\":" << result.ThreadID << ",";
    json << "\n\t\t\t\"ts\":" << result.Start.count();
    json << "\n\t\t}";

    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_current_session)
    {
        m_output_stream << json.str();
        m_output_stream.flush();
    }
}

void ProfilerSingleton::WriteHeader()
{
    m_output_stream << "{\n\t\"otherData\": {},\n\t\"traceEvents\":[\n\t\t{}";
    m_output_stream.flush();
}

void ProfilerSingleton::WriteFooter()
{
    m_output_stream << "\n\t]\n}";
    m_output_stream.flush();
}

void ProfilerSingleton::InternalEndSession()
{
    if (m_current_session)
    {
        WriteFooter();
        m_output_stream.close();
        delete m_current_session;
        m_current_session = nullptr;
    }
}