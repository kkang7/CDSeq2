

#ifndef NECSIM_XOROSHIRO256PLUS_H
#define NECSIM_XOROSHIRO256PLUS_H

void splitMixSeed(uint64_t seed);
uint64_t splimixnext();
void xoshiro256seed(uint64_t seed);
uint64_t xoshiro256rng();

#endif //NECSIM_XOROSHIRO256PLUS_H
