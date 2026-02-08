#pragma once

#include "profiler.h"

#include "func_sig.h"

enum TimeRecordType : uint8_t
{
    None,
    List,
    Accumulate
};

class TimerSingleton
{
public:
    TimerSingleton(const TimerSingleton&) = delete;
    TimerSingleton(TimerSingleton&&)      = delete;

    static TimerSingleton& Get()
    {
        static TimerSingleton instance;
        return instance;
    }

    void Upload(const std::string& name, std::chrono::microseconds time, TimeRecordType record_type, bool output);

    void GetStat(const std::string& name, int last_num, double& avg, double& std);

    double GetAcc(const std::string& name);
    void   clearAcc();

private:
    TimerSingleton()
    {
        if (MPIUtils::is_mpi())
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        }
    }

    int mpi_rank = 0, mpi_size = 1;

    std::unordered_map<std::string, std::vector<std::chrono::microseconds>> m_history_list;
    std::unordered_map<std::string, double>                                 m_history_acc;
};

class Timer
{
public:
    Timer(const std::string& name, TimeRecordType record_type, bool output)
        : m_name(name)
        , m_stopped(false)
        , m_record_type(record_type)
        , m_output(output)
    {
        m_stopped = false;

        size_t pos = m_name.find(cdel_str);
        if (pos != std::string::npos)
        {
            m_name.erase(pos, cdel_str.length());
        }
        m_start_timepoint = std::chrono::steady_clock::now();
    }

    ~Timer()
    {
        if (!m_stopped)
            Stop();
    }

    void Stop()
    {
        auto end_timepoint  = std::chrono::steady_clock::now();
        auto high_res_start = FloatingPointMicroseconds {m_start_timepoint.time_since_epoch()};
        auto elapsed_time =
            std::chrono::time_point_cast<std::chrono::microseconds>(end_timepoint).time_since_epoch() -
            std::chrono::time_point_cast<std::chrono::microseconds>(m_start_timepoint).time_since_epoch();

        TimerSingleton::Get().Upload(m_name, elapsed_time, m_record_type, m_output);
        ProfilerSingleton::Get().Upload({m_name, high_res_start, elapsed_time, std::this_thread::get_id()});

        m_stopped = true;
    }

private:
    std::string                                        m_name;
    std::chrono::time_point<std::chrono::steady_clock> m_start_timepoint;
    bool                                               m_stopped;
    TimeRecordType                                     m_record_type;
    bool                                               m_output = false;
};

#define SCOPE_TIMER_LINE(name, line, type, output)      Timer timer##line(name, type, output)
#define SCOPE_TIMER_PASTELINE(name, line, type, output) SCOPE_TIMER_LINE(name, line, type, output)
#define SCOPE_TIMER(name, type, output)                 SCOPE_TIMER_PASTELINE(name, __LINE__, type, output)
#define FUNCTION_TIMER(type, output)                    SCOPE_TIMER(FUNC_SIG, type, output)
