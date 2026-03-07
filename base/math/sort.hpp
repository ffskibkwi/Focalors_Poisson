#pragma once

template<typename ValueType>
void reorder(ValueType* array, ValueType* buffer, int* indices, int length)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < length; ++i)
    {
        buffer[i] = array[indices[i]];
    }
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < length; ++i)
    {
        array[i] = buffer[i];
    }
}
