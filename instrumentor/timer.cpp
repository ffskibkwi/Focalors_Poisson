#include "timer.h"
#include <cmath>

void TimerSingleton::Upload(const std::string& name, std::chrono::microseconds time)
{
    if (m_regs_stat.find(name) != m_regs_stat.end())
    {
        if (m_history.find(name) != m_history.end())
        {
            m_history[name].push_back(time);
        }
        else
        {
            m_history.insert({name, std::vector<std::chrono::microseconds> {time}});
        }
    }

    if (m_regs_std_cout.find(name) != m_regs_std_cout.end() && mpi_rank == 0)
    {
        if (m_stdcout_enable)
        {
            std::cout << name << " " << static_cast<double>(time.count()) / 1000000.0 << "s" << std::endl;
        }
    }
}

void TimerSingleton::GetStat(const std::string& name, int last_num, double& avg, double& std)
{
    if (m_regs_stat.find(name) == m_regs_stat.end())
    {
        avg = -1.0;
        std = -1.0;

        return;
    }

    std::vector<std::chrono::microseconds>& time_arr = m_history[name];

    // Calculate the average
    double sum = 0.0;

    int begin = std::max(0, (int)time_arr.size() - last_num);
    int count = (int)time_arr.size() - begin;
    for (int i = (int)time_arr.size() - 1; i >= begin; i--)
    {
        sum += static_cast<double>(time_arr[i].count()) / 1000.0;
    }
    avg = sum / count;

    // Calculate the standard deviation
    double variance_sum = 0.0;
    for (int i = (int)time_arr.size() - 1; i >= begin; i--)
    {
        double diff = static_cast<double>(time_arr[i].count()) / 1000.0 - avg;
        variance_sum += diff * diff;
    }
    std = std::sqrt(variance_sum / count);
}
