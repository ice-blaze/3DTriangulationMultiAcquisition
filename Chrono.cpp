#include "Chrono.h"
#include <iostream>
#include <ctime>

using namespace std;

Chrono::Chrono()
{
    restart();

}

Chrono::~Chrono()
{
    //dtor
}

void Chrono::restart()
{
    start = time(NULL);
}
void Chrono::printTime()
{
    int minutes = (getTimeSeconde()/60)+0.5f;
    int secondes = getTimeSeconde()%60;
    cout << "Time : " << minutes<< "m"<<secondes<<"s"<<endl;
//cout << time(NULL) << endl;
}

int Chrono::getTimeSeconde()
{
    return time(NULL)-start;
}

