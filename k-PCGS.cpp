#pragma comment(linker, "/STACK:102400000,102400000")
#include <bits/stdc++.h>

using namespace std;
typedef pair <int, int> P;
typedef pair <double, int> Pd;
const int N = 1e7 + 5;
const int Maxn = 1e7 + 5;
const int Maxk = 105;
const int inf = 0x3f3f3f3f;
map <int, int> anc[N];
vector <int> vec[N], vec2[N];
vector <double> sp_feq_UB[N], sp_feq_LB[N], Feq[N];
int K, dist[N] = {}, out_degree[N], in_degree[N];
bool del[N], del_can[N];
double feq[N], smy_UB[N], smy_LB[N];
int n, m, num, cnt = 0, root = 1;
map <int, int>::iterator ite;
struct Greedy{
    vector <int> vec[Maxn];
    vector <double> sp_feq[Maxn];
    int feq[Maxn], del_can[Maxn], K;
    int n, m, num, cnt, root;
    map <int, int>::iterator ite;
    double add_score(int x, int dist){
        double ans = 0;
        for (int i = 0; i < sp_feq[x].size(); i++)
            ans += 1.0 / (dist + i) * sp_feq[x][i];
        return ans;
    }
    double bfs(int root, int dis[]){
        queue <P> q;
        q.push(P(root, 0));
        bool vis[N] = {};
        vis[root] = true;
        double score = 0;
        while(!q.empty()){
            P p = q.front();
            q.pop();
            int x = p.first, dep = p.second + 1;
            if (dis[x] == inf) score += add_score(x, dep); else score += add_score(x, dep) - add_score(x, dis[x] + 1);
            for (int i = 0; i < vec[x].size(); i++){
                int y = vec[x][i];
                if (vis[y] || dis[y] <= dep) continue;
                q.push(P(y, dep));
                vis[y] = true;
            }
        }
        //cout << root << " " << score << endl;
        return score;
    }
    void update(int root, int dis[]){
        queue <P> q;
        q.push(P(root, 0));
        bool vis[N] = {};
        vis[root] = true;
        dis[root] = 0;
        while(!q.empty()){
            P p = q.front();
            q.pop();
            int x = p.first, dep = p.second + 1;
            //cout << x << " " << dis[x] << endl;
            for (int i = 0; i < vec[x].size(); i++){
                int y = vec[x][i];
                if (vis[y] || dis[y] <= dep) continue;
                dis[y] = dep;
                q.push(P(y, dep));
                vis[y] = true;
            }
        }
    }
    int dis[Maxn], Rank[Maxn], Rank2[Maxn], cntg, cntrank;
    vector <int> solve(int K){
        double ans = 0;
        vector <int> S;
        double score = 0;
        memset(dis, 0x3f, sizeof dis);
        priority_queue <Pd> pq;
        for (int i = 1; i <= n; i++)
            if (!del_can[i]) pq.push(Pd(bfs(i, dis), i));
        while(S.size() < K){
            double max_delta = pq.top().first;
            int x = pq.top().second;
            pq.pop();
            cntg++;
            if (Rank[x] == 0) Rank[x] = ++cntrank;
            double tmp = bfs(x, dis);
            if (max_delta > tmp){
                pq.push(P(tmp, x));
                continue;
            }
            //cout << max_delta << " " << x << endl;
            score += tmp;
            S.push_back(x);
            update(x, dis);
        }
        //sort(S.begin(), S.end());
        cout << cntg << endl;
        cout << score << endl;
        return S;
    }
}Gr;
void get_dist(int root){
    memset(dist, 0, sizeof dist);
    queue <int> q;
    q.push(root), dist[root] = 1;
    while(!q.empty()){
        int x = q.front();
        q.pop();
        for (int i = 0; i < vec[x].size(); i++){
            int y = vec[x][i];
            if (!dist[y]){
                dist[y] = dist[x] + 1;
                q.push(y);
            }
        }
    }
    return;
}
double discount(int x){
    return 1.0 / (x + 1);
}
void del_node1(int K, double UB){
    int du[N] = {}, del_du[N] = {};
    priority_queue <Pd, vector<Pd>, greater<Pd> > pq;
    set <Pd> topk;
    for (int i = 1; i <= n; i++)
        if (topk.size() < K) topk.insert(Pd(1.0 * feq[i] / 2, i));
        else topk.insert(Pd(1.0 * feq[i] / 2, i)), topk.erase(topk.begin());
    for (int i = 1; i <= n; i++){
        du[i] = del_du[i] = out_degree[i];
        del[i] = del_can[i] = 0;
        smy_LB[i] = 1.0 * feq[i] / 2;
        smy_UB[i] = 1.0 * feq[i];
        sp_feq_LB[i].clear();
        sp_feq_LB[i].push_back(feq[i]);
        sp_feq_UB[i].clear();
        sp_feq_UB[i].push_back(feq[i]);
        Feq[i].clear();
        Feq[i].push_back(feq[i]);
        if (du[i] == 0) pq.push(Pd(smy_UB[i], i));
    }
    while(!pq.empty()){
        int x = pq.top().second;
        double val = pq.top().first;
        pq.pop();
        if (val < topk.begin()->first){
            del_can[x] = true;
            if ((vec2[x].size() == 1 || val == 0) && del_du[x] == 0) del[x] = true;
        }
        for (int j = 0; j < vec2[x].size(); j++){
            int y = vec2[x][j];
            du[y]--;
            double tmp_UB = 0, tmp_LB = 0;
            for (int i = 0; i < sp_feq_UB[x].size(); i++)
                if (i + 1 < sp_feq_UB[y].size()) sp_feq_UB[y][i + 1] += sp_feq_UB[x][i];
                else sp_feq_UB[y].push_back(sp_feq_UB[x][i]);
            if (del[x]){
                del_du[y]--;
                for (int i = 0; i < Feq[x].size(); i++)
                    if (i + 1 < Feq[y].size()) Feq[y][i + 1] += Feq[x][i];
                    else Feq[y].push_back(Feq[x][i]);
            }
            if (del_can[x] && vec2[x].size() == 1){
                for (int i = 0; i < sp_feq_LB[x].size(); i++)
                    if (i + 1 < sp_feq_LB[y].size()) sp_feq_LB[y][i + 1] += sp_feq_LB[x][i];
                    else sp_feq_LB[y].push_back(sp_feq_LB[x][i]);
            }
            for (int i = 0; i < sp_feq_UB[y].size(); i++)
                tmp_UB += sp_feq_UB[y][i] * discount(i);
            for (int i = 0; i < sp_feq_LB[y].size(); i++)
                tmp_LB += sp_feq_LB[y][i] * (discount(i) - discount(i + 1));
            if (topk.find(Pd(smy_LB[y], y)) != topk.end()){
                topk.erase(Pd(smy_LB[y], y)), topk.insert(Pd(tmp_LB, y));
            }else{
                topk.insert(Pd(tmp_LB, y)), topk.erase(topk.begin());
            }
            smy_LB[y] = tmp_LB, smy_UB[y] = tmp_UB;
            if (du[y] == 0) pq.push(Pd(tmp_UB, y));
        }
        sp_feq_LB[x].clear();
        sp_feq_UB[x].clear();
        if (del[x]) Feq[x].clear();
    }
    return;
}
void input_greedy(int id[], int re_id[]){
    Gr.n = Gr.m = 0;
    //memset(Gr.del_can, 0, sizeof Gr.del_can);
    for (int i = 1; i <= n; i++)
        if (!del[i]){
            id[i] = ++Gr.n, re_id[Gr.n] = i;
            Gr.del_can[id[i]] = del_can[i];
        }
    for (int i = 1; i <= n; i++){
        if (!del[i]){
            int x = i;
            for (int j = 0; j < vec[x].size(); j++){
                int y = vec[x][j];
                if (!del[y]){
                    //cout << id[x] << " " << id[y] << endl;
                    Gr.m++, Gr.vec[id[x]].push_back(id[y]);
                }
            }
            //cout << id[x] << " " << Feq[x] << endl;
            Gr.sp_feq[id[x]] = Feq[x];
        }
    }
}
double get_score(vector <int> S){
    double score = 0;
    int dist[N] = {};
    queue <int> q;
    for (int i = 0; i < S.size(); i++) q.push(S[i]), dist[S[i]] = 1;
    while(!q.empty()){
        int x = q.front();
        q.pop();
        score += 1.0 * feq[x] / dist[x];
        for (int i = 0; i < vec[x].size(); i++){
            int y = vec[x][i];
            if (!dist[y]){
                dist[y] = dist[x] + 1;
                q.push(y);
            }
        }
    }
    return score;
}
void solve(int K){
    memset(del, 0, sizeof del);
    get_dist(1);
    del_node1(K, 0);
    int cnt1 = 0, cnt2 = 0, id[N] = {}, re_id[N];
    for (int i = 1; i <= n; i++){
        if (!del[i]) cnt1++;
        if (!del_can[i]) cnt2++;
    }
    input_greedy(id, re_id);
    //cout << n << " " << cnt1 << " " << cnt2 << " " << Gr.m << endl;
    vector <int> S = Gr.solve(K);
    for (int i = 0; i < S.size(); i++){
        S[i] = re_id[S[i]];
    }
    cout << get_score(S) << endl;
    //return cnt;
}
void input(){
    //cout << cnt++ << endl;
    //FILE *fin, *fout;
    //fin = fopen("tree-MED.txt", "r");
    //fin = fopen("test_tree.txt", "r");
    //freopen("data.in", "r", stdin);
    //freopen("dp.out", "w", stdout);
    scanf("%d %d", &n, &m);
    for (int i = 1; i <= n; i++) vec[i].clear(), vec2[i].clear(), in_degree[i] = out_degree[i] = 0;
    for (int i = 1; i <= m; i++)
	{
	    int x, y;
		scanf("%d %d", &x, &y);
		vec[x].push_back(y);
		vec2[y].push_back(x);
		out_degree[x]++, in_degree[y]++;
		//fa[mp_id[y]] = mp_id[x];
		//cout << mp_id[x] << " " << mp_id[y] << endl;
	}
	memset(feq, 0, sizeof feq);
	int T;
	scanf("%d", &T);
	//int sum = 0;
	for (int i = 0; i < T; i++)
	{
	    int x, y;
		scanf("%d %d", &x, &y);
		//if (mp_id[x] == 0 && y != 0) cout << "fuck" << endl;
		if (y != 0) {
            feq[x] = y;
            //I.push_back(mp_id[x]);
		}
		//sum += y;
	}
	//sort(I.begin(), I.end(), cmp);
}
int my_main()
{
//    int size = 512 << 20; // 256MB
//    char *p = (char*)malloc(size) + size;
//    __asm__("movl %0, %%esp\n" :: "r"(p));
    //freopen("dag6.txt", "r", stdin);
    //freopen("dag6_vdag.txt", "w", stdout);
    input();
    K = 10;
    //solve(K);
    clock_t t1;
    double t = 0;
    t1 = clock();
    solve(K);
    t = (clock() - t1) * 1.0 / CLOCKS_PER_SEC;
    printf("%.3f\n", t);
    return 0;
}

/*
5 4
1 2
1 3
3 4
3 5
5
1 10
2 20
3 30
4 40
5 30
*/
