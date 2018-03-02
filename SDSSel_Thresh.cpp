#include<bits/stdc++.h>
#include <iostream>
#include <limits.h>
#include <string.h>
#include <queue>
#define V 102
#define REP(x,y,z) for(int x=y;x<=z;x++)
#define FORD(x,y,z) for(int x=y;x>=z;x--)
#define MSET(x,y) memset(x,y,sizeof(x))
#define FOR(x,y) for(__typeof(y.begin()) x=y.begin();x!=y.end();x++)
#define F first
#define S second
#define MP make_pair
#define PB push_back
#define SZ size()
#define M
void RI(){}
template<typename... T>
void RI( int& head, T&... tail ) {
    scanf("%d",&head);
    RI(tail...);
}
using namespace std;
typedef long long LL;
clock_t st,ed;
time_t stwall,edwall;
bool CUT1=false, CUT2=false;
bool SKIP_POST=false;
bool SKIP_DIAMETER=false;
bool SKIP_IMP = false;
int H,P;
void readarg(int argc,char *argv[])
{
	H = P = -1;
	REP(i,0,argc-1)
	{
		if(i<argc-1 && !strcmp(argv[i], "-h"))
		{
			H = atoi(argv[i+1]);
		}
		if(i<argc-1 && !strcmp(argv[i], "-p"))
		{
			P = atoi(argv[i+1]);
		}

		if(!strcmp(argv[i], "-prune1"))
		{
			CUT1 = true;
		}
		if(!strcmp(argv[i], "-prune2"))
		{
			CUT2 = true;
		}
		if(!strcmp(argv[i], "-skip-post"))
		{
			SKIP_POST = true;
		}
		if(!strcmp(argv[i], "-skip-diameter"))
		{
			SKIP_DIAMETER = true;
		}
		if(!strcmp(argv[i], "-skip-imp"))
		{
			SKIP_IMP = true;
		}
	}
}
/**
 * set intersection
 *
 * @param v1,v2: sorted set
 * @return: a sorted array
 */
vector<int> instersection(vector<int> &v1, vector<int> &v2)
{
	vector<int> v3;
	//sort(v1.begin(), v1.end());
	//sort(v2.begin(), v2.end());
	set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
	return v3;
}
int v,c,n;
vector<bool> del,apx,res_cut,ndel;
vector<vector<pair<int,double>>> e;
vector<vector<double>> ediv;
vector<pair<int,int>> deg;
vector<int> dis;
bool iter = false;
/**
 * initializa and read graph input
 */
void init()
{
	int x,y;
	int e1,e2,e3;
	double z;

	scanf("%d %d %d %d %d",&v,&c,&e1,&e2,&e3);
	n = v+c;
	e = vector<vector<pair<int,double>>> (n+1);
	ediv = vector<vector<double>> (c+1);
	REP(i,0,c) ediv[i] = vector<double>(c+1);
	deg = vector<pair<int,int>> (n+1);
	del = vector<bool> (n+1);
	apx = vector<bool> (n+1);
  res_cut = vector<bool> (V-2);
  ndel = vector<bool> (V-2);

	REP(i,1,e1) //read social edges
	{
		scanf("%d %d",&x,&y);
		e[x].PB( MP(y,0.0) );
		e[y].PB( MP(x,0.0) );
	}
	REP(i,1,e2) //read diveristy edgess
	{
		scanf("%d %d %lf",&x,&y,&z);
		e[x].PB( MP(y,z) );
		e[y].PB( MP(x,z) );
		//ediv[x].PB( MP(y,z) );
		//ediv[y].PB( MP(x,z) );
		//printf("%d %d\n",x-v,y-v);
		ediv[x-v][y-v] = ediv[y-v][x-v] = z;
	}
	REP(i,1,e3) //read channel edge
	{
		scanf("%d %d",&x,&y);
		e[x].PB( MP(y,0.0) );
		e[y].PB( MP(x,0.0) );
	}
}
/**
 * get degree of all vertices(excepted deleted nodes)
 * deg[i].first: number of viewer connected to i
 * deg[i].second: number of channel conntected to i
 */
void get_degree()
{
	REP(i,1,n) deg[i] = MP(0,0);
	REP(i,1,n)if(!del[i]) for(auto j:e[i])if(!del[j.F]) //foreach edge i->j
	{
		if(j.F>v) deg[i].S++;
		else deg[i].F++;
	}
}
/**
 * return "cur" 's h-hop set
 *
 * @param cur: a vertex in viewers set
 * @return pair<vector<int>,vector<int>>: first vector are viewers in h-hop set,
 * second vector are channels in h-hop set(channels doesn't satisfy the constraint are deleted)
 */
pair<vector<int>,vector<int>> get_h2(int cur)
{
	queue<int> q;
  vector<int> r1,r2;
	vector<int> dis2(n+1);
	dis2[cur]=1;//distance from cur

	//find social vertices
	for(auto i:e[cur]) if(i.F<=v && !dis2[i.F] && !del[i.F])
	{
		dis2[i.F] = dis2[cur]+1;
		q.push(i.F);
	}
	while(!q.empty())//BFS queue
	{
		cur = q.front();
		q.pop();
		if(dis2[cur]>=H+1) continue;

		for(auto i:e[cur]) if(i.F<=v && !dis2[i.F] && !del[i.F])
		{
			dis2[i.F] = dis2[cur]+1;
			q.push(i.F);
		}
	}

	//channels
  REP(i,1,v) if(dis2[i]>0 && dis2[i]<=H+1) r1.PB(i);
	//channels
  if( (int)r1.size()<P ) return MP(r1,r2);
	vector<int> cnt(n+1);
	REP(i,v+1,v+c) for(auto j:e[i]) if(j.F<=v && dis2[j.F]) cnt[i]++;
	REP(i,v+1,v+c) if(cnt[i]>=P && !del[i]) r2.PB(i);
	return MP(r1,r2);
}
//given channels，return their preference users
vector<int> get_viewers(vector<int> &ch)//
{
	vector<int> re;
	vector<bool> viss(n+1);
	for(int i:ch) for(auto j:e[i]) if(j.F<=v && !del[j.F]) viss[j.F]=true;
	REP(i,1,v) if(viss[i]) re.PB(i);
	return re;
}
//given viewers，find their preference channels (doesn't delete channel < p preference)
vector<int> get_channels(vector<int> &vi)//
{
	vector<int> re;
	vector<bool> viss(n+1);
	for(int i:vi) for(auto j:e[i]) if(j.F>v && !del[j.F]) viss[j.F]=true;
	REP(i,v+1,v+c) if(viss[i]) re.PB(i);
	return re;
}
//given channels，find the one with max total diversity(and must be in set y)
int max_channel_diversity(vector<int> &x,vector<int> &y)//
{
	int re=-1;
	double mx = -1e9;

	vector<double> sum(c+1);
	int sz = (int)x.size();
	REP(i,0,sz-1)REP(j,i+1,sz-1) if(x[i]>v && x[j]>v)
	{
		sum[x[i]-v] += ediv[x[i]-v][x[j]-v];
		sum[x[j]-v] += ediv[x[i]-v][x[j]-v];
	}

	for(int i:y) if(sum[i-v]>mx)
	{
		mx = sum[i-v];
		re = i;
	}
	return re;
}
//given channels，find the total diversity of this subgraph
double total_diversity(vector<int> &x) //
{
	double re=0;
	for(int i:x) if(!del[i] && i>v) for(int j:x) if(!del[j] && j>v) re += ediv[i-v][j-v];
	return re;
}
//given one channel, find total diversity
double vertex_total_diversity(int vv) //
{
	double re=0;
	int i=vv;
	REP(j,v+1,v+c) if(!del[j]) re += ediv[i-v][j-v];
	return re;
}
//given channel set, find objective value
double get_objective(vector<int> &x)//
{
	if((int)x.size()==0) return 0.0;

	double re=0.0;
	for(int i:x) if(!del[i] && i>v) for(int j:x) if(!del[j] && j>v) re += ediv[i-v][j-v];
	return re / (double)x.size();
}
//find shortest path from node "root"(on viewers graph，nodes deleted are also considered)
void shortest_path(int root)
{
	int INF = 100000000;
	int cur;
	queue<int> q;

    //BFS
	dis = vector<int> (n+1, INF);
	{
		dis[root] = 0;
		q.push(root);
		while(!q.empty())
		{
			cur = q.front();
			q.pop();

			for(auto j:e[cur])if(j.F<=v && dis[j.F]==INF)
			{
				dis[j.F] = dis[cur]+1;
				q.push(j.F);
			}
		}
	}
}

/* Returns true if there is a path from source 's' to sink 't' in
  residual graph. Also fills parent[] to store the path */
int bfs(double rGraph[V][V], int s, int t, int parent[])
{
    // Create a visited array and mark all vertices as not visited
    bool visited[V];
    memset(visited, 0, sizeof(visited));

    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    // Standard BFS Loop
    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (int v=0; v<V; v++)
        {
            if (visited[v]==false && rGraph[u][v] > 0)
            {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == true);
}

// A DFS based function to find all reachable vertices from s.  The function
// marks visited[i] as true if i is reachable from s.  The initial values in
// visited[] must be false. We can also use BFS to find reachable vertices
void dfs(double rGraph[V][V], int s, bool visited[])
{
    visited[s] = true;
    for (int i = 0; i < V; i++)
       if (rGraph[s][i]>0 && !visited[i])
           dfs(rGraph, i, visited);
}

// Prints the minimum s-t cut
void minCut(double graph[V][V], int s, int t)
{
    int u, v;

    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
    double rGraph[V][V]; // rGraph[i][j] indicates residual capacity of edge i-j
    for (u = 0; u < V; u++)
        for (v = 0; v < V; v++)
             rGraph[u][v] = graph[u][v];

    int parent[V];  // This array is filled by BFS and to store path

    // Augment the flow while tere is path from source to sink
    while (bfs(rGraph, s, t, parent))
    {
        // Find minimum residual capacity of the edhes along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        int path_flow = INT_MAX;
        for (v=t; v!=s; v=parent[v])
        {
            u = parent[v];
            path_flow = min(double(path_flow), rGraph[u][v]);
        }

        // update residual capacities of the edges and reverse edges
        // along the path
        for (v=t; v != s; v=parent[v])
        {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
    }

    // Flow is maximum now, find vertices reachable from s
    bool visited[V];
    memset(visited, false, sizeof(visited));
    dfs(rGraph, s, visited);
    for (int i = 0; i < V-2; i++)
    {
      if(!visited[i]&&!ndel[i])
      {
        iter = true;
        ndel[i] = true;
      }

    }
    //**** delete all the nodes from res_cut = true
    // Print all edges that are from a reachable vertex to
    // non-reachable vertex in the original graph
    //for (int i = 0; i < V; i++)
      //for (int j = 0; j < V; j++)
        // if (visited[i] && !visited[j] && graph[i][j]>0)
          //    cout << i << " - " << j << endl;

    //return;
}

void sdssel()
{
    //delete channel < p viewers
	get_degree();
	REP(i,v+1,v+c) if(deg[i].F<P) del[i]=true;

	//delete viewer with no preference edge
	get_degree();
	REP(i,1,v) if(deg[i].S<=0) del[i]=true;

	double kapxobj = 0.0;
	double kapximpobj = 0.0;
	int vb=-1, qb=-1;
	vector<int> h2v,h2c;
	vector<int> kapx, fapx;
	vector<int> kapximp, fapximp;
	vector<int> cfit;

	REP(Vs,1,v) if(!del[Vs]) //select v belong to I, I = I - {v}
	{
		/*//new cut added
		vector<int> vvtmp2; vvtmp2.PB(Vs);
		vector<int> cctmp2 = get_channels(vvtmp2);//C(v)
		//q_bar = argmax_{q belong to C(v)} (Delta C)
		qb = max_channel_diversity(cctmp2, cctmp2);
		if( vertex_total_diversity(qb) <= 2*kapxobj )
		{
			continue;
		}*/

        //get H_v^2(V) and H_v^2(C)

		tie(h2v,h2c) = get_h2(Vs);
		if( (int)h2v.size()<P || (int)h2c.size()<2) continue;
    iter = true;
    double graph[V][V];
    double lambda = 0.0;
    vector<double> indiv(V-2);
    REP(i,0,V-3)
    {
      ndel[i] = true;
      for(int ip:h2c) if (ip==i+1) ndel[i] = false;
    }
    while(iter)
    {
      iter = false;
      REP(i,0,V-3)
      {
        REP(j,0,V-3)
        {
          if(!ndel[i] && !ndel[j])
          {
            graph[i][j] = ediv[i+1-v][j+1-v];
            graph[j][i] = ediv[i+1-v][j+1-v];
            indiv[i] = indiv[i]+graph[i][j];
            lambda = lambda + graph[i][j];
          }
          else
          {
            graph[i][j] = 0;
            graph[j][i] = 0;
          }
        }
      }
      int k=0;
      REP(i,0,V-3) if (!ndel[i]) k=k+1;
      if(k>0)lambda = lambda / (2*k);
          //prune 2
      REP(i,0,V-3)
      {
        if (indiv[i]>=lambda && !ndel[i])
        {
          graph[i][V-2]=indiv[i]-lambda;
          graph[V-2][i]=indiv[i]-lambda;
          graph[i][V-1]=0;
          graph[V-1][i]=0;
        }
        if (indiv[i]<lambda && !ndel[i])
        {
          graph[i][V-1]=lambda-indiv[i];
          graph[V-1][i]=lambda-indiv[i];
          graph[i][V-2]=0;
          graph[V-2][i]=0;
        }
        if (ndel[i])
        {
          graph[i][V-1]=0;
          graph[V-1][i]=0;
          graph[i][V-2]=0;
          graph[V-2][i]=0;
        }
      }
      graph[V-1][V-2]=0;
      graph[V-2][V-1]=0;
      minCut(graph, V-2, V-1);

    }
    vector<int> newh2c;
    REP(i,0,V-3) if (!ndel[i]) newh2c.PB(i+1);
    h2c = newh2c;
		if(CUT2)
		{
			vector<int> ptmp; ptmp.PB(Vs);
			vector<int> vtmp = get_channels(ptmp);
			if((int)(instersection(h2c, vtmp).size())==0)
				continue;
		}
		//prune 1
		if(CUT1 && vb!=-1)
		{
			vector<int> ptmp; ptmp.PB(vb);
			vector<int> vtmp = get_channels(ptmp);

			if( total_diversity(h2c) <= 2*get_objective(vtmp) )
				continue;
		}

		//calc diversity
		vector<double> h2cdiv = vector<double>(n+1);//h2c div = diversity of each node in H_v^2(C)
		for(int i:h2c) for(int j:h2c) if(i>v && j>v) h2cdiv[i] += ediv[i-v][j-v];

		double h2ctotdiv = total_diversity(h2c);
		double h2cobj = h2ctotdiv / (int)h2c.size();
    if( kapx.size()==0 || h2cobj > kapxobj )
    {
      vector<int> viewer = get_viewers(h2c);
      vector<int> tmpf = instersection(h2v, viewer);

      if((int)tmpf.size() >= P)
      {
        kapx = h2c;
        kapxobj = h2cobj;
        fapx = instersection(h2v, viewer);
        vb = Vs;
      }
    }
		//while(h2c.size())
		//{
			//if( kapx.size()==0 || h2cobj > kapxobj )
			//{
				//vector<int> viewer = get_viewers(h2c);
				//vector<int> tmpf = instersection(h2v, viewer);

				//if((int)tmpf.size() >= P)
				//{
					//kapx = h2c;
					//kapxobj = h2cobj;
					//fapx = instersection(h2v, viewer);
					//vb = Vs;
				//}
			//}
			//q_bar = min_channel_diversity(h2c)
			//qb = -1;
			//for(int i:h2c) if(qb==-1 || h2cdiv[i]<h2cdiv[qb]) qb = i;

			//delete node q_bar
			////maintain diversity after q_bar deleted
			//for(int i:h2c) h2cdiv[i] -= ediv[i-v][qb-v];

			//maintain h2c's total diversity
			//for(int i:h2c) h2ctotdiv -= 2*ediv[i-v][qb-v];

			//delete q_bar from h2c
			//h2c.erase(remove(h2c.begin(), h2c.end(), qb), h2c.end());

			//recalc h2c's objective value
			//h2cobj = h2ctotdiv / (int)h2c.size();
		//}
	}

	//IMP
	if(!SKIP_IMP)
	{
    tie(h2v,h2c) = get_h2(vb);
    for(int i:h2c) if(lower_bound(kapx.begin(), kapx.end(), i) == kapx.end()) cfit.PB(i);
		kapximp = kapx;
		kapximpobj = kapxobj;
		fapximp = fapx;

		//delete from cfit
		vector<int> tp;
		for(int i:cfit)
		{
			double s=0;
			for(auto j:e[i]) if(!del[j.F] && j.F>v && lower_bound(h2c.begin(), h2c.end(), j.F)!=h2c.end())
				s += j.S;
			if(s < kapxobj)
				tp.PB(i);
		}
		for(int i:tp) cfit.erase(remove(cfit.begin(), cfit.end(), i), cfit.end());

        //record diversity of every nodes in kapximp
        vector<double> kapximpdiv(n+1);
		for(int i:kapximp) for(int j:kapximp) if(i>v && j>v) kapximpdiv[i] += ediv[i-v][j-v];
		while(kapximpobj<=kapxobj && cfit.size())
		{
			//qb = max_channel_diversity(kapximp, cfit);
			qb = -1;
			for(int i:cfit) if(qb==-1 || kapximpdiv[i]>kapximpdiv[qb]) qb = i;
			vector<int> tt2; tt2.PB(qb);
			vector<int> viewer = get_viewers( tt2 );

			//if qb not in kapximp, insert it
			if( lower_bound(kapximp.begin(),kapximp.end(),qb) == kapximp.end() )
			{
				kapximp.PB(qb);
				sort(kapximp.begin(), kapximp.end());
				kapximp.resize( unique(kapximp.begin(), kapximp.end()) - kapximp.begin() );
				//recalc objective value
				kapximpobj = get_objective(kapximp);

				kapximpdiv = vector<double>(n+1);
				for(int i:kapximp) for(int j:kapximp) if(i>v && j>v) kapximpdiv[i] += ediv[i-v][j-v];
			}

			//fapximp
			for(int i:viewer) if(lower_bound(h2v.begin(), h2v.end(), i) != h2v.end()) fapximp.PB(i);
			sort(fapximp.begin(), fapximp.end());
			fapximp.resize( unique(fapximp.begin(), fapximp.end()) - fapximp.begin() );

			cfit.erase(remove(cfit.begin(), cfit.end(), qb), cfit.end());

			//if(get_objective(kapximp)>get_objective(kapx))
			if(kapximpobj > kapxobj)
			{
				kapx = kapximp;
				kapxobj = kapximpobj;
				fapx = fapximp;
			}
		}
	}

	//POST process
	vector<bool> delf(n+1); //delete from apx
	if(!SKIP_POST)
	{
		for(int i:fapx) if(!delf[i])
		{
		    //check shortest path
			shortest_path(i);
			bool flag=false;
			for(int j:fapx) if(!delf[j]) if(dis[j]>H) { flag=true; break; } //distance > H
			if(!flag) continue;

			//reduce channel
			vector<int> cnt(n+1);
			vector<int> newch;
			vector<bool> vnow(n+1);
			for(int j:fapx) if(!delf[j] && j!=i) vnow[j]=true;

			//count the viewers
			REP(j,v+1,v+c) for(auto k:e[j]) if(vnow[k.F]) cnt[j]++;
			for(int j:kapx) if(cnt[j]>=P) newch.PB(j);

			double newchobj = get_objective(newch);
			//if(get_objective(newch) - get_objective(kapx) >= -1e-7)
			if(newchobj - kapxobj >= -1e-7) //if objective value won't become smaller, delete the node
			{
				delf[i] = true;
				kapx = newch;
				kapxobj = newchobj;
			}
		}
	}

	vector<int> newf;
	for(int i:fapx) if(!delf[i]) newf.PB(i);
	fapx = newf;

	//degree
	fill(del.begin(), del.end(), true);
	for(int i:fapx) del[i] = false;
	for(int i:kapx) del[i] = false;
	get_degree();

	//DELETE viewer without preference
	REP(i,1,v) if(!del[i] && deg[i].S<1)
		delf[i] = true;
	newf.clear();
	for(int i:fapx) if(!delf[i]) newf.PB(i);
	fapx = newf;

	fill(del.begin(), del.end(), true);
	for(int i:fapx) del[i] = false;
	for(int i:kapx) del[i] = false;
	get_degree();

	ed = clock();//end timer
  edwall = time(NULL);

	//OUTPUT
	//calc farest pair
	int farest = 0;
	if(!SKIP_DIAMETER)
	{
		for(int i:fapx)
		{
			shortest_path(i);
			for(int j:fapx) farest = max(farest, dis[j]);
		}
	}

    //check is fesible solution(satisfy degree constraint)
	bool fes = true;
	REP(i,1,v) if(!del[i] && deg[i].S<1) fes = false;
	REP(i,v+1,v+c) 	if(!del[i] && deg[i].F<P) fes = false;

    //average viewer of each channel
	double avg = 0.0;
	REP(i,v+1,v+c) if(!del[i]) avg += deg[i].F;
	avg /= (double)kapx.size();

	if((int)kapx.size()==0)
	{
		fes = false;
		avg = 0.0;
	}
	printf("CPU Time %f, Wall clock time %f, Viewer V cap %d, %d %d %f %d %f\n",
		(double)(ed-st)/CLOCKS_PER_SEC,
    difftime(edwall, stwall),
    (int)vb,
		(int)(fes && farest<=H),
		(int)fes,
		get_objective(kapx),
		farest,
		avg
	);
	printf("%d %d\n", (int)fapx.size(), (int)kapx.size());
	for(int i:fapx) printf("%d ",i); puts("");
	for(int i:kapx) printf("%d ",i); puts("");
}
int main(int argc,char *argv[])
{
	st = clock();
  stwall = time(NULL);
	readarg(argc,argv);
	init();
	sdssel();
	return 0;
}




// Driver program to test above functions
//int main()
//{
    // Let us create a graph shown in the above example
    //int graph[V][V] = { {0, 16, 13, 0, 0, 0},
                        //{0, 0, 10, 12, 0, 0},
                        //{0, 4, 0, 0, 14, 0},
                        //{0, 0, 9, 0, 0, 20},
                        //{0, 0, 0, 7, 0, 4},
                        //{0, 0, 0, 0, 0, 0}
                      //};

    //minCut(graph, 0, 5);

    //return 0;
//}
