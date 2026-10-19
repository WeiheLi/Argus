#include "Pandora.hpp"
#include <math.h>
extern unsigned long g_hash_seed_offset;
extern int g_pandora_def;

Pandora::Pandora(int depth, int width, int lgn)
{

	pandora_.depth = depth;
	pandora_.width = width;
	pandora_.lgn = lgn;
	pandora_.sum = 0;
	pandora_.counts = new SBucket *[depth * width];
	for (int i = 0; i < depth * width; i++)
	{
		pandora_.counts[i] = (SBucket *)calloc(1, sizeof(SBucket));
		memset(pandora_.counts[i], 0, sizeof(SBucket));
		pandora_.counts[i]->key[0] = '\0';
	}

	pandora_.hash = new unsigned long[depth];
	pandora_.scale = new unsigned long[depth];
	pandora_.hardner = new unsigned long[depth];
	char name[] = "Pandora";
	unsigned long seed = g_hash_seed_offset + AwareHash((unsigned char *)name, strlen(name), 13091204281, 228204732751, 6620830889);
	for (int i = 0; i < depth; i++)
	{
		pandora_.hash[i] = GenHashSeed(seed++);
	}
	for (int i = 0; i < depth; i++)
	{
		pandora_.scale[i] = GenHashSeed(seed++);
	}
	for (int i = 0; i < depth; i++)
	{
		pandora_.hardner[i] = GenHashSeed(seed++);
	}
}

Pandora::~Pandora()
{
	for (int i = 0; i < pandora_.depth * pandora_.width; i++)
	{
		free(pandora_.counts[i]);
	}
	delete[] pandora_.hash;
	delete[] pandora_.scale;
	delete[] pandora_.hardner;
	delete[] pandora_.counts;
}

void Pandora::Update(unsigned char *key, val_tp val)
{
	unsigned long bucket = 0;
	int keylen = pandora_.lgn / 8;
	pandora_.sum += 1;
	Pandora::SBucket *sbucket;
	int flag = 0;
	int value = 0;
	int default_value = g_pandora_def;
	long min = 99999999;
	int loc = -1;
	int k;
	int index;
	for (int i = 0; i < pandora_.depth; i++)
	{
		bucket = MurmurHash64A(key, keylen, pandora_.hardner[i]) % pandora_.width;
		index = i * pandora_.width + bucket;
		sbucket = pandora_.counts[index];
		if (sbucket->key[0] == '\0' && sbucket->count == 0 && sbucket->status == 0)
		{
			memcpy(sbucket->key, key, keylen);
			flag = 1;
			sbucket->status = 1;
			sbucket->count = 1;
			return;
		}
		else if (memcmp(key, sbucket->key, keylen) == 0)
		{
			if (sbucket->status == 0)
			{
				flag = 1;
				sbucket->count += 1;
				sbucket->status = 1;
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
		sbucket = pandora_.counts[loc];
		if (sbucket->status == 1)
		{
			return;
		}

		if (default_value - sbucket->coldcount < 1)
			value = 1;
		else
			value = default_value - sbucket->coldcount;

		k = rand() % (int)((sbucket->count * value) + 1.0) + 1.0;
		if (k > (int)((sbucket->count * value)))
		{
			memcpy(sbucket->key, key, keylen);
			sbucket->count = 1;
			sbucket->coldcount = 0;
			sbucket->status = 1;
		}
	}
}

void Pandora::Query(val_tp thresh, std::vector<std::pair<key_tp, val_tp>> &results)
{
	for (int i = 0; i < pandora_.width * pandora_.depth; i++)
	{

		if (pandora_.counts[i]->count > (int)thresh)
		{
			key_tp reskey;
			memcpy(reskey.key, pandora_.counts[i]->key, pandora_.lgn / 8);
			std::pair<key_tp, val_tp> node;
			node.first = reskey;
			node.second = pandora_.counts[i]->count;
			results.push_back(node);
		}
	}
}

void Pandora::Inactivity()
{
	for (int i = 0; i < pandora_.depth * pandora_.width; i++)
	{
		if (pandora_.counts[i]->status == 0)
		{
			pandora_.counts[i]->coldcount = pandora_.counts[i]->coldcount + 1;
		}
		if (pandora_.counts[i]->status == 1)
		{
			pandora_.counts[i]->coldcount = pandora_.counts[i]->coldcount - 1;
			if (pandora_.counts[i]->coldcount < 0)
			{
				pandora_.counts[i]->coldcount = 0;
			}
		}
	}
}

void Pandora::NewWindow()
{
	for (int i = 0; i < pandora_.depth * pandora_.width; i++)
	{
		pandora_.counts[i]->status = 0;
	}
}

val_tp Pandora::PointQuery(unsigned char *key)
{
	return Low_estimate(key);
}

val_tp Pandora::Low_estimate(unsigned char *key)
{

	val_tp ret = 0, max = 0, min = 999999999;
	for (int i = 0; i < pandora_.depth; i++)
	{
		unsigned long bucket = MurmurHash64A(key, pandora_.lgn / 8, pandora_.hardner[i]) % pandora_.width;

		unsigned long index = i * pandora_.width + bucket;
		if (memcmp(pandora_.counts[index]->key, key, pandora_.lgn / 8) == 0)
		{
			max += pandora_.counts[index]->count;
		}
		index = i * pandora_.width + (bucket + 1) % pandora_.width;
		if (memcmp(key, pandora_.counts[i]->key, pandora_.lgn / 8) == 0)
		{
			max += pandora_.counts[index]->count;
		}
	}
	return max;
}

val_tp Pandora::Up_estimate(unsigned char *key)
{

	val_tp ret = 0, max = 0, min = 999999999;
	for (int i = 0; i < pandora_.depth; i++)
	{
		unsigned long bucket = MurmurHash64A(key, pandora_.lgn / 8, pandora_.hardner[i]) % pandora_.width;

		unsigned long index = i * pandora_.width + bucket;
		if (memcmp(pandora_.counts[index]->key, key, pandora_.lgn / 8) == 0)
		{
			max += pandora_.counts[index]->count;
		}
		if (pandora_.counts[index]->count < min)
			min = pandora_.counts[index]->count;
		index = i * pandora_.width + (bucket + 1) % pandora_.width;
		if (memcmp(key, pandora_.counts[i]->key, pandora_.lgn / 8) == 0)
		{
			max += pandora_.counts[index]->count;
		}
	}
	if (max)
		return max;
	return min;
}

val_tp Pandora::GetCount()
{
	return pandora_.sum;
}

void Pandora::Reset()
{
	pandora_.sum = 0;
	for (int i = 0; i < pandora_.depth * pandora_.width; i++)
	{
		pandora_.counts[i]->count = 0;
		memset(pandora_.counts[i]->key, 0, pandora_.lgn / 8);
	}
}

void Pandora::SetBucket(int row, int column, val_tp sum, long count, unsigned char *key)
{
	int index = row * pandora_.width + column;
	pandora_.counts[index]->count = count;
	memcpy(pandora_.counts[index]->key, key, pandora_.lgn / 8);
}

Pandora::SBucket **Pandora::GetTable()
{
	return pandora_.counts;
}
