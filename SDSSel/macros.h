#ifndef MACROS_H
#define MACROS_H
#define REP(x, y, z) for (int x = y; x <= z; x++)
#define FORD(x, y, z) for (int x = y; x >= z; x--)
#define MSET(x, y) memset(x, y, sizeof(x))
#define FOR(x, y) for (__typeof(y.begin()) x = y.begin(); x != y.end(); x++)
#define F first
#define S second
#define MP std::make_pair
#define PB push_back
#define SZ size()
#define M
#endif
