#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
using namespace std;
#define NUM 5000
#define SEED 3

int main(int argc, char *argv[])
{
    ofstream MyFile("baza5000.txt");
    //int tc=strtol(argv[1],NULL,10);
    srand(SEED);
    //#pragma omp parallel for num_threads(tc)
    for(int i=0;i<NUM;i++)
    {
        for(int j=0;j<NUM;j++)
        {
            MyFile<<((double)((rand()%100)*pow(-1.0,(double)(rand()%2))))/10;
            MyFile<<" ";
        }
        MyFile<<endl;
    }
    for(int i =0;i<NUM;i++)
    {
        MyFile<<(((double)(rand()%200))/10)*(-1);
        MyFile<<" ";
    }
    MyFile<<endl;
    for(int i =0;i<NUM;i++)
    {
        MyFile<<(200+rand()%9800);
        MyFile<<" ";
    }

    return 0;
}
