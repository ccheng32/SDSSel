#include <bits/stdc++.h>
#include <stdlib.h>
#include <unistd.h>
#include <mutex>
#include <thread>
#define REP(x, y, z) for (int x = y; x <= z; x++)
#define FORD(x, y, z) for (int x = y; x >= z; x--)
#define MSET(x, y) memset(x, y, sizeof(x))
#define FOR(x, y) for (__typeof(y.begin()) x = y.begin(); x != y.end(); x++)
#define F first
#define S second
#define MP make_pair
#define PB push_back
#define SZ size()
#define M
void RI() {}
template <typename... T>
void RI(int &head, T &... tail) {
  scanf("%d", &head);
  RI(tail...);
}
using namespace std;
typedef long long LL;
clock_t st, ed;
time_t stwall, edwall;
struct timespec stthread, edthread;
std::mutex muapx;
bool CUT1 = false, CUT2 = false;
bool SKIP_POST = false;
bool SKIP_DIAMETER = false;
bool SKIP_IMP = false;
int H, P, Cr, nbrV;
double timespec2msec(const struct timespec *ts) {
  return (double)ts->tv_sec + (double)ts->tv_nsec * 1e-9;
}
void readarg(int argc, char *argv[]) {
  H = P = -1;
  REP(i, 0, argc - 1) {
    if (i < argc - 1 && !strcmp(argv[i], "-h")) {
      H = atoi(argv[i + 1]);
    }
    if (i < argc - 1 && !strcmp(argv[i], "-p")) {
      P = atoi(argv[i + 1]);
    }
    if (i < argc - 1 && !strcmp(argv[i], "-c")) {
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
}
/**
 * 集合交集
 *
 * @param v1,v2: sorted set
 * @return: a sorted array
 */
vector<int> instersection(vector<int> &v1, vector<int> &v2) {
  vector<int> v3;
  // sort(v1.begin(), v1.end());
  // sort(v2.begin(), v2.end());
  set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                   back_inserter(v3));
  return v3;
}
int v, c, n;
vector<bool> del, apx;
vector<vector<pair<int, double>>> e;
vector<vector<double>> ediv;
vector<pair<int, int>> deg;
vector<int> dis;
vector<int> Proc_threads;
vector<int> sharedkapx;
int sharedkapxobj;
vector<int> sharedfapx;
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
  e = vector<vector<pair<int, double>>>(n + 1);
  ediv = vector<vector<double>>(c + 1);
  REP(i, 0, c) ediv[i] = vector<double>(c + 1);
  deg = vector<pair<int, int>>(n + 1);
  del = vector<bool>(n + 1);
  apx = vector<bool>(n + 1);

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
 * @return pair<vector<int>,vector<int>>: 第一個vector是h-hop裡面的viewer,
 * 第二個是對應到的channel(有刪掉不符合constraint的)
 */
pair<vector<int>, vector<int>> get_h2(int cur) {
  queue<int> q;
  vector<int> r1, r2;
  vector<int> dis2(n + 1);
  dis2[cur] = 1;

  // find social vertices
  for (auto i : e[cur])
    if (i.F <= v && !dis2[i.F] && !del[i.F]) {
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
  vector<int> cnt(n + 1);
  REP(i, v + 1, v + c) for (auto j : e[i]) if (j.F <= v && dis2[j.F]) cnt[i]++;
  REP(i, v + 1, v + c) if (cnt[i] >= P && !del[i]) r2.PB(i);
  return MP(r1, r2);
}
//傳channel進去，找出他們的preference users
vector<int> get_viewers(vector<int> &ch)  //
{
  vector<int> re;
  vector<bool> viss(n + 1);
  for (int i : ch)
    for (auto j : e[i])
      if (j.F <= v && !del[j.F]) viss[j.F] = true;
  REP(i, 1, v) if (viss[i]) re.PB(i);
  return re;
}
//傳一些viewers進去，找出他們有的channels (沒有刪除小於p個人訂閱的)
vector<int> get_channels(vector<int> &vi)  //
{
  vector<int> re;
  vector<bool> viss(n + 1);
  for (int i : vi)
    for (auto j : e[i])
      if (j.F > v && !del[j.F]) viss[j.F] = true;
  REP(i, v + 1, v + c) if (viss[i]) re.PB(i);
  return re;
}
//傳一些channel進去，找出裡面total diversity最大的(必須在y set裡面)
int max_channel_diversity(vector<int> &x, vector<int> &y)  //
{
  int re = -1;
  double mx = -1e9;

  vector<double> sum(c + 1);
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
double total_diversity(vector<int> &x)  //
{
  double re = 0;
  //	vector<bool> visx(n+1);
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
double get_objective(vector<int> &x)  //
{
  if ((int)x.size() == 0) return 0.0;

  double re = 0.0;
  //	vector<bool> visx(n+1);
  //	vector<double> sum(n+1);
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
  queue<int> q;

  dis = vector<int>(n + 1, INF);
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
void shared_apx(vector<int> kapxthread, vector<int> fapxthread,
                int kapxthreadobj, int vbthread, int Proc, int num) {
  std::lock_guard<std::mutex> guard(muapx);
  Proc_threads[num] = Proc;
  if (kapxthreadobj > sharedkapxobj) {
    sharedkapx = kapxthread;
    sharedkapxobj = kapxthreadobj;
    sharedfapx = fapxthread;
    sharedvb = vbthread;
  }
}
void sdsselthread(int startrange, int endrange, int threadnum) {
  double kapxobj = 0.0;
  int vb = -1, qb = -1;
  vector<int> h2v, h2c, h2v1, h2c1;
  vector<int> kapx, fapx;
  int Proc = 0;
  // vector<int> kapximp, fapximp;
  // vector<int> cfit;

  REP(Vs, startrange, endrange) if (!del[Vs]) {
    /*//new cut added
    vector<int> vvtmp2; vvtmp2.PB(Vs);
    vector<int> cctmp2 = get_channels(vvtmp2);
    qb = max_channel_diversity(cctmp2, cctmp2);
    if( vertex_total_diversity(qb) <= 2*kapxobj )
    {
            continue;
    }*/

    tie(h2v, h2c) = get_h2(Vs);

    if ((int)h2v.size() < P) continue;
    tie(h2v1, h2c1) = get_h2(Vs);
    Proc = Proc + 1;
    if (CUT2) {
      vector<int> ptmp;
      ptmp.PB(Vs);
      vector<int> vtmp = get_channels(ptmp);
      if ((int)(instersection(h2c, vtmp).size()) == 0) continue;
    }
    if (CUT1 && vb != -1) {
      vector<int> ptmp;
      ptmp.PB(vb);
      vector<int> vtmp = get_channels(ptmp);

      if (total_diversity(h2c) <= 2 * get_objective(vtmp)) continue;
    }

    // calc diversity
    vector<double> h2cdiv = vector<double>(n + 1);
    for (int i : h2c)
      for (int j : h2c)
        if (i > v && j > v) h2cdiv[i] += ediv[i - v][j - v];

    double h2ctotdiv = total_diversity(h2c);
    double h2cobj = h2ctotdiv / (int)h2c.size();
    // switch kapx and h2c
    while (h2c.size()) {
      if (kapx.size() == 0 || h2cobj > kapxobj) {
        vector<int> viewer = get_viewers(h2c);
        vector<int> tmpf = instersection(h2v, viewer);

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
  int j;
  // REP(j,0,endrange) j++;
  shared_apx(kapx, fapx, kapxobj, vb, Proc, threadnum);
}
// IMP
// sort(h2c.begin(), h2c.end());
// sort(kapx.begin(), kapx.end());
void IMP() {
  double kapximpobj = 0.0;
  vector<int> kapximp, fapximp;
  vector<int> cfit;
  vector<int> kapx = sharedkapx;
  int kapxobj = sharedkapxobj;
  vector<int> fapx = sharedfapx;
  int vb = sharedvb;
  vector<int> h2v, h2c;
  tie(h2v, h2c) = get_h2(vb);
  if (!SKIP_IMP) {
    for (int i : h2c)
      if (lower_bound(kapx.begin(), kapx.end(), i) == kapx.end()) cfit.PB(i);
    kapximp = kapx;
    kapximpobj = kapxobj;
    fapximp = fapx;

    // delete from cfit
    vector<int> tp;
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

    vector<double> kapximpdiv(n + 1);
    for (int i : kapximp)
      for (int j : kapximp)
        if (i > v && j > v) kapximpdiv[i] += ediv[i - v][j - v];
    while (kapximpobj <= kapxobj && cfit.size()) {
      // qb = max_channel_diversity(kapximp, cfit);
      int qb = -1;
      for (int i : cfit)
        if (qb == -1 || kapximpdiv[i] > kapximpdiv[qb]) qb = i;
      vector<int> tt2;
      tt2.PB(qb);
      vector<int> viewer = get_viewers(tt2);

      // kapximp
      if (lower_bound(kapximp.begin(), kapximp.end(), qb) == kapximp.end()) {
        kapximp.PB(qb);
        sort(kapximp.begin(), kapximp.end());
        kapximp.resize(unique(kapximp.begin(), kapximp.end()) -
                       kapximp.begin());
        kapximpobj = get_objective(kapximp);

        kapximpdiv = vector<double>(n + 1);
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
  vector<bool> delf(n + 1);  // delete from apx
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
        vector<int> cnt(n + 1);
        vector<int> newch;
        vector<bool> vnow(n + 1);
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

  vector<int> newf;
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
      for (int j : fapx) farest = max(farest, dis[j]);
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
int main(int argc, char *argv[]) {  // Now you have to multithread this part
  st = clock();
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stthread);
  stwall = time(NULL);
  readarg(argc, argv);
  init();
  preprocessing();
  nbrV = 0;
  Proc_threads = vector<int>(Cr);
  // int nbrthread=std::thread::hardware_concurrency();
  REP(i, 1, v) if (!del[i]) nbrV++;
  int nbrVThread = (int)nbrV / Cr;
  int remain = nbrV % Cr;
  vector<thread> threads(Cr - 1);

  int startrange = 1;
  int endrange = 1;
  for (int i = 0; i < Cr - 1; ++i) {
    int k = 0;
    int j = startrange;
    int limitthread;
    if (remain > 0) {
      limitthread = nbrVThread + 1;
      remain--;
    } else {
      limitthread = nbrVThread;
    }
    while (k < limitthread) {
      if (!del[j]) k++;
      j++;
    }
    endrange = j;
    threads[i] = thread(sdsselthread, startrange, endrange, i);
    startrange = endrange + 1;
  }
  sdsselthread(startrange, v, Cr - 1);
  std::for_each(threads.begin(), threads.end(),
                std::mem_fn(&std::thread::join));
  IMP();
  return 0;
}
