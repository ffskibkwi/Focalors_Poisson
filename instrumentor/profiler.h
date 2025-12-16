#pragma once

#include "poisson_base/base/parallel/mpi/mpi_misc.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

    extern const std::string cdel_str;

using FloatingPointMicroseconds = std::chrono::duration<double, std::micro>;

struct ProfileResult
{
    std::string               Name;
    FloatingPointMicroseconds Start;
    std::chrono::microseconds ElapsedTime;
    std::thread::id           ThreadID;
};

struct ProfilerSession
{
    std::string Name;
};

class ProfilerSingleton
{
public:
    ProfilerSingleton(const ProfilerSingleton&) = delete;
    ProfilerSingleton(ProfilerSingleton&&)      = delete;

    static ProfilerSingleton& Get()
    {
        static ProfilerSingleton instance;
        return instance;
    }

    void BeginSession(const std::string& name, const std::string& filepath = "results.json");

    void EndSession();

    void Upload(const ProfileResult& result);

    void EnableWritingProfile() { is_writing_enabled = true; }

    void DisableWritingProfile() { is_writing_enabled = false; }

    void RegisterNameToCalcAverageFlops(const std::string& name)
    {
        m_registered_names_to_calc_average_flops.insert(name);
    }

private:
    ProfilerSingleton()
        : m_current_session(nullptr)
    {
        if (MPIUtils::is_mpi())
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        }
    }

    ~ProfilerSingleton() { EndSession(); }

    void WriteHeader();

    void WriteFooter();

    void InternalEndSession();

private:
    int mpi_rank = 0, mpi_size = 1;

    bool             is_writing_enabled = false;
    std::mutex       m_mutex;
    ProfilerSession* m_current_session;
    std::ofstream    m_output_stream;

    std::unordered_set<std::string> m_registered_names_to_calc_average_flops;
};

#define PROFILE_BEGIN_SESSION(name, filepath) ProfilerSingleton::Get().BeginSession(name, filepath)
#define PROFILE_END_SESSION()                 ProfilerSingleton::Get().EndSession()
