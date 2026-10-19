#include "PSketch.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;
extern int g_psketch_k;

PSketch::PSketch(int depth, int width, int lgn) {

	psketch_.depth = depth;
	psketch_.width = width;
	psketch_.lgn = lgn;
	psketch_.sum = 0;
	psketch_.counts = new SBucket *[depth*width];
	for (int i = 0; i < depth*width; i++) {
		psketch_.counts[i] = (SBucket*)calloc(1, sizeof(SBucket));
		memset(psketch_.counts[i], 0, sizeof(SBucket));
		psketch_.counts[i]->key[0] = '\0';
	}

	psketch_.hash = new unsigned long[depth];
	psketch_.scale = new unsigned long[depth];
	psketch_.hardner = new unsigned long[depth];
	char name[] = "PSketch";
	unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char*)name, strlen(name), 13091204281, 228204732751, 6620830889);
	for (int i = 0; i < depth; i++) {
		psketch_.hash[i] = GenHashSeed(seed++);
	}
	for (int i = 0; i < depth; i++) {
		psketch_.scale[i] = GenHashSeed(seed++);
	}
	for (int i = 0; i < depth; i++) {
		psketch_.hardner[i] = GenHashSeed(seed++);
	}
}

PSketch::~PSketch() {
	for (int i = 0; i < psketch_.depth*psketch_.width; i++) {
		free(psketch_.counts[i]);
	}
	delete[] psketch_.hash;
	delete[] psketch_.scale;
	delete[] psketch_.hardner;
	delete[] psketch_.counts;
}

void PSketch::Update(unsigned char* key, val_tp val) {
	unsigned long bucket = 0;
	int keylen = psketch_.lgn / 8;psketch_.sum += 1;
	PSketch::SBucket *sbucket;
	int flag = 0;
	long min = 99999999; int loc = -1; int k; int index;
	for (int i = 0; i < psketch_.depth; i++) {
		bucket = MurmurHash64A(key, keylen, psketch_.hardner[i]) % psketch_.width;
		index = i * psketch_.width + bucket;
		sbucket = psketch_.counts[index];
		if (sbucket->key[0] == '\0'&&sbucket->count==0&&sbucket->status==0) { 
			memcpy(sbucket->key, key, keylen);
			flag = 1;
			sbucket->status = 1;
			sbucket->count = 1;
                       return;
		}
		else if (memcmp(key, sbucket->key, keylen) == 0) {
                       if (sbucket->status==0)
			{
                          flag = 1;
			   sbucket->count += 1;
                          sbucket->status=1; 
                        }
                       return;
		}
        if (sbucket->count < min)
		{
			min = sbucket->count;
			loc = index; 
		}
	}
	if (flag == 0 && loc >= 0) 
	{
		sbucket = psketch_.counts[loc]; 
                if(sbucket->status==1) {return;}
		k = rand() % (int)(g_psketch_k*(sbucket->count+sbucket->hotcount)+1.0) + 1.0;
		if (k > (int)(g_psketch_k*(sbucket->count+sbucket->hotcount))) 
		{                    
                       memcpy(sbucket->key,key,keylen);
			sbucket->count += 1;
			sbucket->status=1;
		}
	}
}

void PSketch::Query(val_tp thresh, std::vector<std::pair<key_tp, val_tp> >&results) {
	myset res; 
	for (int i = 0; i < psketch_.width*psketch_.depth; i++) {
                     
		if (psketch_.counts[i]->count > (int)thresh) {
			key_tp reskey;
			memcpy(reskey.key, psketch_.counts[i]->key, psketch_.lgn / 8);
			res.insert(reskey);

		}

	}
	for (auto it = res.begin(); it != res.end(); it++) {
		val_tp resval = 99999999; val_tp max = 0;
		for (int j = 0; j < psketch_.depth; j++) {
			unsigned long bucket = MurmurHash64A((*it).key, psketch_.lgn / 8, psketch_.hardner[j]) % psketch_.width;
			unsigned long index = j * psketch_.width + bucket;
			if (memcmp(psketch_.counts[index]->key, (*it).key, psketch_.lgn / 8) == 0) {
				max += psketch_.counts[index]->count;			
			}
			
		}
		if (max != 0)resval = max;

			key_tp key;
			memcpy(key.key, (*it).key, psketch_.lgn / 8);
			std::pair<key_tp, val_tp> node;
			node.first = key;
			node.second = resval;
			results.push_back(node);
	}
}


void PSketch::QueryHotness(){ 
for (int i = 0; i < psketch_.depth*psketch_.width; i++){
                        if(psketch_.counts[i]->status == 0)
                         {
                           psketch_.counts[i]->hotcount = psketch_.counts[i]->hotcount - 1;
                           if(psketch_.counts[i]->hotcount<0)
                            {psketch_.counts[i]->hotcount = 0;}
                         }
                        if(psketch_.counts[i]->status == 1)
                          {psketch_.counts[i]->hotcount = psketch_.counts[i]->hotcount + 1;}
	}
}

void PSketch::NewWindow(){ 
for (int i = 0; i < psketch_.depth*psketch_.width; i++){
			psketch_.counts[i]->status = 0;
	}
}

val_tp PSketch::PointQuery(unsigned char* key) {
	return Low_estimate(key);
}

val_tp PSketch::Low_estimate(unsigned char* key) {


val_tp ret = 0, max = 0, min = 999999999;
	for (int i = 0; i < psketch_.depth; i++) {
		unsigned long bucket = MurmurHash64A(key, psketch_.lgn / 8, psketch_.hardner[i]) % psketch_.width;

		unsigned long index = i * psketch_.width + bucket;
		if (memcmp(psketch_.counts[index]->key, key, psketch_.lgn / 8) == 0)
		{
			max += psketch_.counts[index]->count;
		}
		index = i * psketch_.width + (bucket + 1) % psketch_.width;
		if (memcmp(key, psketch_.counts[i]->key, psketch_.lgn / 8) == 0)
		{
			max += psketch_.counts[index]->count;
		}

	}
	return max;


}

val_tp PSketch::Up_estimate(unsigned char* key) {

val_tp ret = 0, max = 0, min = 999999999;
	for (int i = 0; i < psketch_.depth; i++) {
		unsigned long bucket = MurmurHash64A(key, psketch_.lgn / 8, psketch_.hardner[i]) % psketch_.width;

		unsigned long index = i * psketch_.width + bucket;
		if (memcmp(psketch_.counts[index]->key, key, psketch_.lgn / 8) == 0)
		{
			max += psketch_.counts[index]->count;
		}
		if (psketch_.counts[index]->count < min)min = psketch_.counts[index]->count;
		index = i * psketch_.width + (bucket + 1) % psketch_.width;
		if (memcmp(key, psketch_.counts[i]->key, psketch_.lgn / 8) == 0)
		{
			max += psketch_.counts[index]->count;
		}

	}
	if (max)return max;
	return min;

}

val_tp PSketch::GetCount() {
	return psketch_.sum;
}

void PSketch::Reset() {
	psketch_.sum = 0;
	for (int i = 0; i < psketch_.depth*psketch_.width; i++) {
		psketch_.counts[i]->count = 0;
		memset(psketch_.counts[i]->key, 0, psketch_.lgn / 8);
	}
}

void PSketch::SetBucket(int row, int column, val_tp sum, long count, unsigned char* key) {
	int index = row * psketch_.width + column;
	psketch_.counts[index]->count = count;
	memcpy(psketch_.counts[index]->key, key, psketch_.lgn / 8);
}

PSketch::SBucket** PSketch::GetTable() {
	return psketch_.counts;
}






