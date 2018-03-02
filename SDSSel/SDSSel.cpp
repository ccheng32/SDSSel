#include <bits/stdc++.h>
#include <stdlib.h>
#include <unistd.h>
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
void RI() {}
template <typename... T>
void RI(int &head, T &... tail) {
  scanf("%d", &head);
  RI(tail...);
}

typedef long long LL;
clock_t st, ed;
time_t stwall, edwall;
struct timespec stthread, edthread;
bool CUT1 = false, CUT2 = false;
bool SKIP_POST = false;
bool SKIP_DIAMETER = false;
bool SKIP_IMP = false;

// std::maximum number of hops between any two users
int H = -1;

// minimum number of channels each user is interested in
int P = -1;

// number of threads used
int Cr = 1;

double timespec2msec(const struct timespec *ts) {
  return (double)ts->tv_sec + (double)ts->tv_nsec * 1e-9;
}

#define READ_ARG_SUCCESS 0xABCD
#define READ_ARG_FAIL 0xBCDA
int readarg(int argc, char *argv[]) {
  H = P = -1;
  REP(i, 0, argc - 1) {
    if (!strcmp(argv[i], "-h")) {
      H = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "-p")) {
      P = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "-c")) {
      Cr = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "-prune1")) {
      CUT1 = true;
    }
    if (!strcmp(argv[i], "-prune2")) {
      CUT2 = true;
    }
    if (!strcmp(argv[i], "-skip-post")) {
      SKIP_POST = true;
    }
    if (!strcmp(argv[i], "-skip-diameter")) {
      SKIP_DIAMETER = true;
    }
    if (!strcmp(argv[i], "-skip-imp")) {
      SKIP_IMP = true;
    }
  }

  if (H <= 0 || P <= 0) {
    printf("Please provide parameters for hops (-h) and pref (-p)\n");
    return READ_ARG_FAIL;
  } else {
    return READ_ARG_SUCCESS;
  }
}

/**
 * 集合交集
 *
 * @param v1,v2: sorted set
 * @return: a sorted array
 */
std::vector<int> instersection(std::vector<int> &v1, std::vector<int> &v2) {
  std::vector<int> v3;
  // sort(v1.begin(), v1.end());
  // sort(v2.begin(), v2.end());
  set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                   back_inserter(v3));
  return v3;
}

int v, c, n;
std::vector<bool> del, apx;
std::vector<std::vector<std::pair<int, double>>> e;
std::vector<std::vector<double>> ediv;
std::vector<std::pair<int, int>> deg;
std::vector<int> dis;
std::vector<int> sharedkapx;
int sharedkapxobj;
std::vector<int> sharedfapx;
int sharedvb;
/**
 * 初始化、讀圖檔
 */
void init() {
  int x, y;
  int e1, e2, e3;
  double z;

  scanf("%d %d %d %d %d", &v, &c, &e1, &e2, &e3);
  n = v + c;
  e = std::vector<std::vector<std::pair<int, double>>>(n + 1);
  ediv = std::vector<std::vector<double>>(c + 1);
  REP(i, 0, c) ediv[i] = std::vector<double>(c + 1);
  deg = std::vector<std::pair<int, int>>(n + 1);
  del = std::vector<bool>(n + 1);
  apx = std::vector<bool>(n + 1);

  REP(i, 1, e1) {
    scanf("%d %d", &x, &y);
    e[x].PB(MP(y, 0.0));
    e[y].PB(MP(x, 0.0));
  }
  REP(i, 1, e2) {
    scanf("%d %d %lf", &x, &y, &z);
    e[x].PB(MP(y, z));
    e[y].PB(MP(x, z));
    // ediv[x].PB( MP(y,z) );
    // ediv[y].PB( MP(x,z) );
    // printf("%d %d\n",x-v,y-v);
    ediv[x - v][y - v] = ediv[y - v][x - v] = z;
  }
  REP(i, 1, e3) {
    scanf("%d %d", &x, &y);
    e[x].PB(MP(y, 0.0));
    e[y].PB(MP(x, 0.0));
  }
}

/**
 * 算所有點(除了被deleted的點)的degree
 * deg[i].first: i 連到的 viewer 數量
 * deg[i].second: i 連到的 channel 數量
 */
void get_degree() {
  REP(i, 1, n) deg[i] = MP(0, 0);
  REP(i, 1, n) if (!del[i]) for (auto j : e[i]) if (!del[j.F]) {
    if (j.F > v)
      deg[i].S++;
    else
      deg[i].F++;
  }
}
/**
 * 回傳 cur 的 h hop set
 *
 * @param cur 中心的viewer點
 * @eturn std::pair<std::vector<int>,std::vector<int>>:
 * 第一個std::vector是h-hop裡面的viewer,
 * 第二個是對應到的channel(有刪掉不符合constraint的)
 */
std::pair<std::vector<int>, std::vector<int>> get_h2(int cur) {
  std::queue<int> q;
  std::vector<int> r1, r2;

  // distance from cur to all viewers
  std::vector<int> dis2(n + 1);
  dis2[cur] = 1;

  // get all viewers at most h hops away from cur
  for (auto i : e[cur])
    if (i.F <= v && !del[i.F]) {
      dis2[i.F] = dis2[cur] + 1;
      q.push(i.F);
    }
  while (!q.empty()) {
    cur = q.front();
    q.pop();
    if (dis2[cur] >= H + 1) continue;

    for (auto i : e[cur])
      if (i.F <= v && !dis2[i.F] && !del[i.F]) {
        dis2[i.F] = dis2[cur] + 1;
        q.push(i.F);
      }
  }
  REP(i, 1, v) if (dis2[i] > 0 && dis2[i] <= H + 1) r1.PB(i);
  // channels
  // if( (int)r1.size()<P ) return MP(r1,r2);
  std::vector<int> cnt(n + 1);
  REP(i, v + 1, v + c) for (auto j : e[i]) if (j.F <= v && dis2[j.F]) cnt[i]++;
  REP(i, v + 1, v + c) if (cnt[i] >= P && !del[i]) r2.PB(i);
  return MP(r1, r2);
}
//傳channel進去，找出他們的preference users
std::vector<int> get_viewers(std::vector<int> &ch)  //
{
  std::vector<int> re;
  std::vector<bool> viss(n + 1);
  for (int i : ch)
    for (auto j : e[i])
      if (j.F <= v && !del[j.F]) viss[j.F] = true;
  REP(i, 1, v) if (viss[i]) re.PB(i);
  return re;
}
//傳一些viewers進去，找出他們有的channels (沒有刪除小於p個人訂閱的)
std::vector<int> get_channels(std::vector<int> &vi)  //
{
  std::vector<int> re;
  for (int i : vi)
    for (auto j : e[i])
      if (j.F > v && !del[j.F]) re.PB(j.F);
  return re;
}
//傳一些channel進去，找出裡面total diversity最大的(必須在y set裡面)
int max_channel_diversity(std::vector<int> &x, std::vector<int> &y)  //
{
  int re = -1;
  double mx = -1e9;

  std::vector<double> sum(c + 1);
  int sz = (int)x.size();
  REP(i, 0, sz - 1) REP(j, i + 1, sz - 1) if (x[i] > v && x[j] > v) {
    sum[x[i] - v] += ediv[x[i] - v][x[j] - v];
    sum[x[j] - v] += ediv[x[i] - v][x[j] - v];
  }

  for (int i : y)
    if (sum[i - v] > mx) {
      mx = sum[i - v];
      re = i;
    }
  return re;
}
//傳一些channel進去，找出整張子圖的total diversity
double total_diversity(std::vector<int> &x)  //
{
  double re = 0;
  //	std::vector<bool> visx(n+1);
  //	for(int i:x) visx[i]=true;
  for (int i : x)
    if (!del[i] && i > v)
      for (int j : x)
        if (!del[j] && j > v) re += ediv[i - v][j - v];
  //	for(int i:x) for(auto j:ediv[i]) if(!del[j.F] && visx[j.F] && j.F>v) re
  //+= j.S;

  return re;
}
//傳一個channel進去，找出在整張圖的total diversity
double vertex_total_diversity(int vv)  //
{
  double re = 0;
  int i = vv;  // for(auto j:ediv[i]) if(!del[j.F] && j.F>v) re += j.S;
  REP(j, v + 1, v + c) if (!del[j]) re += ediv[i - v][j - v];
  return re;
}
//找出一個channel set的obj value
double get_objective(std::vector<int> &x)  //
{
  if ((int)x.size() == 0) return 0.0;

  double re = 0.0;
  //	std::vector<bool> visx(n+1);
  //	std::vector<double> sum(n+1);
  //	for(int i:x) visx[i]=true;
  for (int i : x)
    if (!del[i] && i > v)
      for (int j : x)
        if (!del[j] && j > v) re += ediv[i - v][j - v];
  //	for(int i:x) for(auto j:ediv[i]) if(!del[j.F] && visx[j.F] && j.F>v)
  // sum[i] += j.S;
  //	REP(i,v+1,v+c)
  //	{
  //		re += sum[i];
  //	}
  return re / (double)x.size();
}
//求以一個點為出發點的最短路(在viewers原圖上做，會考慮被刪掉的點)
void shortest_path(int root) {
  int INF = 100000000;
  int cur;
  std::queue<int> q;

  dis = std::vector<int>(n + 1, INF);
  {
    dis[root] = 0;
    q.push(root);
    while (!q.empty()) {
      cur = q.front();
      q.pop();

      for (auto j : e[cur])
        if (j.F <= v && dis[j.F] == INF) {
          dis[j.F] = dis[cur] + 1;
          q.push(j.F);
        }
    }
  }
}
void preprocessing() {
  get_degree();
  REP(i, v + 1, v + c) if (deg[i].F < P) del[i] = true;
  get_degree();
  REP(i, 1, v) if (deg[i].S <= 0) del[i] = true;
}
void shared_apx(std::vector<int> kapxthread, std::vector<int> fapxthread,
                int kapxthreadobj, int vbthread) {
  if (kapxthreadobj > sharedkapxobj) {
    sharedkapx = kapxthread;
    sharedkapxobj = kapxthreadobj;
    sharedfapx = fapxthread;
    sharedvb = vbthread;
  }
}
void sdssel() {
  double kapxobj = 0.0;
  int vb = -1, qb = -1;
  std::vector<int> kapx, fapx;

  for (int Vs = 1; Vs <= v; Vs++) {
    std::vector<int> h2v, h2c;
    if (!del[Vs]) {
      std::tie(h2v, h2c) = get_h2(Vs);

      // viewer group is smaller than P, skipped
      if ((int)h2v.size() < P) continue;

      if (CUT2) {
        std::vector<int> ptmp;
        ptmp.PB(Vs);
        std::vector<int> vtmp = get_channels(ptmp);
        if ((int)(instersection(h2c, vtmp).size()) == 0) continue;
      }

      // Channel-set Difference Pruning (CDP)
      if (CUT1 && vb != -1) {
        std::vector<int> ptmp;
        ptmp.PB(vb);
        std::vector<int> vtmp = get_channels(ptmp);

        if (total_diversity(h2c) <= 2 * get_objective(vtmp)) continue;
      }

      // calc diversity
      std::vector<double> h2cdiv = std::vector<double>(n + 1);
      for (int i : h2c)
        for (int j : h2c)
          if (i > v && j > v) h2cdiv[i] += ediv[i - v][j - v];

      double h2ctotdiv = total_diversity(h2c);
      double h2cobj = h2ctotdiv / (int)h2c.size();
      // switch kapx and h2c
      while (h2c.size()) {
        if (kapx.size() == 0 || h2cobj > kapxobj) {
          std::vector<int> viewer = get_viewers(h2c);
          std::vector<int> tmpf = instersection(h2v, viewer);

          if ((int)tmpf.size() >= P) {
            kapx = h2c;
            kapxobj = h2cobj;
            fapx = instersection(h2v, viewer);
            vb = Vs;
          }
        }
        // qb = min_channel_diversity(h2c);
        qb = -1;
        for (int i : h2c)
          if (qb == -1 || h2cdiv[i] < h2cdiv[qb]) qb = i;
        for (int i : h2c) h2cdiv[i] -= ediv[i - v][qb - v];
        for (int i : h2c) h2ctotdiv -= 2 * ediv[i - v][qb - v];
        h2c.erase(remove(h2c.begin(), h2c.end(), qb), h2c.end());
        h2cobj = h2ctotdiv / (int)h2c.size();
      }
    }
  }

  shared_apx(kapx, fapx, kapxobj, vb);
}
// IMP
// sort(h2c.begin(), h2c.end());
// sort(kapx.begin(), kapx.end());
void IMP() {
  double kapximpobj = 0.0;
  std::vector<int> kapximp, fapximp;
  std::vector<int> cfit;
  std::vector<int> kapx = sharedkapx;
  int kapxobj = sharedkapxobj;
  std::vector<int> fapx = sharedfapx;
  int vb = sharedvb;
  std::vector<int> h2v, h2c;
  std::tie(h2v, h2c) = get_h2(vb);
  if (!SKIP_IMP) {
    for (int i : h2c)
      if (lower_bound(kapx.begin(), kapx.end(), i) == kapx.end()) cfit.PB(i);
    kapximp = kapx;
    kapximpobj = kapxobj;
    fapximp = fapx;

    // delete from cfit
    std::vector<int> tp;
    for (int i : cfit) {
      double s = 0;
      for (auto j : e[i])
        if (!del[j.F] && j.F > v &&
            lower_bound(h2c.begin(), h2c.end(), j.F) != h2c.end())
          s += j.S;
      if (s < kapxobj) tp.PB(i);
    }
    for (int i : tp)
      cfit.erase(remove(cfit.begin(), cfit.end(), i), cfit.end());

    std::vector<double> kapximpdiv(n + 1);
    for (int i : kapximp)
      for (int j : kapximp)
        if (i > v && j > v) kapximpdiv[i] += ediv[i - v][j - v];
    while (kapximpobj <= kapxobj && cfit.size()) {
      // qb = std::max_channel_diversity(kapximp, cfit);
      int qb = -1;
      for (int i : cfit)
        if (qb == -1 || kapximpdiv[i] > kapximpdiv[qb]) qb = i;
      std::vector<int> tt2;
      tt2.PB(qb);
      std::vector<int> viewer = get_viewers(tt2);

      // kapximp
      if (lower_bound(kapximp.begin(), kapximp.end(), qb) == kapximp.end()) {
        kapximp.PB(qb);
        sort(kapximp.begin(), kapximp.end());
        kapximp.resize(unique(kapximp.begin(), kapximp.end()) -
                       kapximp.begin());
        kapximpobj = get_objective(kapximp);

        kapximpdiv = std::vector<double>(n + 1);
        for (int i : kapximp)
          for (int j : kapximp)
            if (i > v && j > v) kapximpdiv[i] += ediv[i - v][j - v];
      }

      // fapximp
      for (int i : viewer)
        if (lower_bound(h2v.begin(), h2v.end(), i) != h2v.end()) fapximp.PB(i);
      sort(fapximp.begin(), fapximp.end());
      fapximp.resize(unique(fapximp.begin(), fapximp.end()) - fapximp.begin());

      cfit.erase(remove(cfit.begin(), cfit.end(), qb), cfit.end());

      // if(get_objective(kapximp)>get_objective(kapx))
      if (kapximpobj > kapxobj) {
        kapx = kapximp;
        kapxobj = kapximpobj;
        fapx = fapximp;
      }
    }
  }

  // POST
  std::vector<bool> delf(n + 1);  // delete from apx
  if (!SKIP_POST) {
    for (int i : fapx)
      if (!delf[i]) {
        shortest_path(i);
        bool flag = false;
        for (int j : fapx)
          if (!delf[j])
            if (dis[j] > H) {
              flag = true;
              break;
            }
        if (!flag) continue;

        // reduce channel
        std::vector<int> cnt(n + 1);
        std::vector<int> newch;
        std::vector<bool> vnow(n + 1);
        for (int j : fapx)
          if (!delf[j] && j != i) vnow[j] = true;

        // count the viewers
        REP(j, v + 1, v + c) for (auto k : e[j]) if (vnow[k.F]) cnt[j]++;
        for (int j : kapx)
          if (cnt[j] >= P) newch.PB(j);

        double newchobj = get_objective(newch);
        // if(get_objective(newch) - get_objective(kapx) >= -1e-7)
        if (newchobj - kapxobj >= -1e-7) {
          delf[i] = true;
          kapx = newch;
          kapxobj = newchobj;
        }
      }
  }

  std::vector<int> newf;
  for (int i : fapx)
    if (!delf[i]) newf.PB(i);
  fapx = newf;

  // degree
  fill(del.begin(), del.end(), true);
  for (int i : fapx) del[i] = false;
  for (int i : kapx) del[i] = false;
  get_degree();

  // DELETE viewer without preference
  REP(i, 1, v) if (!del[i] && deg[i].S < 1) delf[i] = true;
  newf.clear();
  for (int i : fapx)
    if (!delf[i]) newf.PB(i);
  fapx = newf;

  fill(del.begin(), del.end(), true);
  for (int i : fapx) del[i] = false;
  for (int i : kapx) del[i] = false;
  get_degree();
  ed = clock();  // end timer
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &edthread);
  edwall = time(NULL);
  // printf("end time %f\n", difftime(edwall, stwall));

  // OUTPUT
  int farest = 0;
  if (!SKIP_DIAMETER) {
    for (int i : fapx) {
      shortest_path(i);
      for (int j : fapx) farest = std::max(farest, dis[j]);
    }
  }

  bool fes = true;
  REP(i, 1, v) if (!del[i] && deg[i].S < 1) fes = false;
  REP(i, v + 1, v + c) if (!del[i] && deg[i].F < P) fes = false;

  double avg = 0.0;
  REP(i, v + 1, v + c) if (!del[i]) avg += deg[i].F;
  avg /= (double)kapx.size();

  if ((int)kapx.size() == 0) {
    fes = false;
    avg = 0.0;
  }
  // printf("%f\n", timespec2msec(&edthread) - timespec2msec(&stthread));
  printf(
      " CPUTime main thread %f, CPUTime %f, Wall clock time %f, Viewer V cap "
      "%d, %d %d %f %d %f\n",
      timespec2msec(&edthread) - timespec2msec(&stthread),
      (double)(ed - st) / CLOCKS_PER_SEC, difftime(edwall, stwall), (int)vb,
      (int)(fes && farest <= H), (int)fes, get_objective(kapx), farest, avg);

  printf("\n------------RESULT-----------------\n");
  printf("number of users:%d\n", (int)fapx.size());
  for (int i : fapx) printf("%d ", i);
  printf("\n");
  printf("number of channels:%d\n", (int)kapx.size());
  for (int i : kapx) printf("%d ", i);
  printf("\n");
}

/*
 * Function: main function
 * -----------------------
 * program entry point
 */
int main(int argc, char *argv[]) {  // Now you have to multithread this part
  st = clock();
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stthread);
  stwall = time(NULL);
  if (readarg(argc, argv) == READ_ARG_FAIL) {
    return 1;
  }

  init();

  preprocessing();

  sdssel();

  IMP();
  return 0;
}
