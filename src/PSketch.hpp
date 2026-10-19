#ifndef PSketch_H
#define PSketch_H
#include <vector>
#include <unordered_set>
#include <utility>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "datatypes.hpp"
extern "C"
{
#include "hash.h"
#include "util.h"
}



class PSketch {

    typedef struct SBUCKET_type {

        short int count;
        short int hotcount;
        unsigned char key[LGN];
        uint8_t status;

    } SBucket;



    struct psketch_type {

        //Counter to count total degree
        val_tp sum;
        //Counter table
        SBucket **counts;

        //Outer sketch depth and width
        int depth;
        int width;

        //# key word bits
        int lgn;

        unsigned long *hash, *scale, *hardner;
    };

    public:
	PSketch(int depth, int width, int lgn);

    ~PSketch();

    void Update(unsigned char* key, val_tp value);

    val_tp PointQuery(unsigned char* key);

    void Query(val_tp thresh, myvector& results);

    void NewWindow();

    void QueryHotness();

    val_tp Low_estimate(unsigned char* key);

    val_tp Up_estimate(unsigned char* key);

    val_tp GetCount();

    void Reset();

    void MergeAll(PSketch** psketch_arr, int size);

    private:

    void SetBucket(int row, int column, val_tp sum, long count, unsigned char* key);

	PSketch::SBucket** GetTable();

	psketch_type psketch_;
};

#endif
