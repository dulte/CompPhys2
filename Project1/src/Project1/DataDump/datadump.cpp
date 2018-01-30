#include "datadump.h"



/*
template<class T>
DataDump<T>::DataDump(std::string m_location, std::string stamp_location)
    :DataDump(m_location)
{
    stampfile->open(stamp_location,std::ofstream::out);
    include_stamp = true;
}*/



template<class T>
void DataDump<T>::push_back(T data_point){
    data.push_back(data_point);
}

template<class T>
void DataDump<T>::dump(T data_point){
    outfile << data_point;
}

template<class T>
void DataDump<T>::dump(T data_point,double stamp){
    dump(data_point);
    if(include_stamp){
        stampfile << stamp;
    }

}

template<class T>
void DataDump<T>::dump_all(){
    int data_size = data.size();
    for(int i = 0; i<data_size;i++){
        outfile << data[i] << "\n";
    }
    if(include_stamp){
        int stamp_size = data_stamp.size();
        for(int i = 0; i<   stamp_size;i++){
            stampfile << data_stamp[i];
        }
    }
}

template class DataDump<double>;
template class DataDump<int>;
template class DataDump<std::string>;
//template class DataDump<std::vector<int>>;
//template class DataDump<std::vector<double>>;





