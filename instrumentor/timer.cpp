#include "timer.h"
#include <cmath>

void TimerSingleton::Upload(const std::string&        name,
                            std::chrono::microseconds time,
                            TimeRecordType            record_type,
                            bool                      output)
{
    double time_s = static_cast<double>(time.count()) / 1000000.0;

    if (record_type == TimeRecordType::List)
    {
        if (m_history_list.find(name) != m_history_list.end())
        {
            m_history_list[name].push_back(time);
        }
        else
        {
            m_history_list.insert({name, std::vector<std::chrono::microseconds> {time}});
        }
    }
    else if (record_type == TimeRecordType::Accumulate)
    {
        if (m_history_acc.find(name) != m_history_acc.end())
        {
            m_history_acc[name] += time_s;
        }
        else
        {
            m_history_acc.insert({name, time_s});
        }
    }

    if (output && mpi_rank == 0)
    {
        std::cout << name << " " << time_s << "s" << std::endl;
    }
}

void TimerSingleton::GetStat(const std::string& name, int last_num, double& avg, double& std)
{
    if (m_history_list.find(name) == m_history_list.end())
    {
        avg = -1.0;
        std = -1.0;

        return;
    }

    std::vector<std::chrono::microseconds>& time_arr = m_history_list[name];

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

double TimerSingleton::GetAcc(const std::string& name)
{
    if (m_history_acc.find(name) != m_history_acc.end())
        return m_history_acc[name];
    else
        return 0.0;
}

void TimerSingleton::clearAcc() { m_history_acc.clear(); }