#include <bits/stdc++.h>
using namespace std;

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
//datatype snippets
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
//stl snippets
typedef pair<int,int> pii;
typedef pair<ll,ll> pll;
typedef vector<bool> vb;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef vector<ld> vld;
typedef vector<string> vs;
typedef vector<pii> vpii;
typedef vector<pll> vpll;
typedef vector<vi> vvi;
typedef vector<vb> vvb;
typedef vector<vector<ll>> vvll;
typedef vector<vpii> vvpii;
typedef vector<vpll> vvpll;
typedef vector<vld> vvld;
typedef set<int> si;
typedef set<ll> sll;
typedef map<int,int> mii;
typedef map<ll,ll> mll;
typedef map<int,pii> mipii;
typedef map<ll,pll> mlpll;
typedef queue<int> qi;
typedef queue<ll> qll;
typedef deque<int> dqi;
typedef deque<ll> dqll;

#define all(v) v.begin(), v.end()
#define eb emplace_back
#define mp make_pair
#define lb lower_bound
#define ub upper_bound
#define fi first
#define se second
//io snippets
#define yes cout<<"YES\n"
#define no cout<<"NO\n"
#define fastIO ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL)
template<typename... Args>
void input(Args&... args){
    ((cin>>args),...);
}
template<typename... Args>
void print(Args... args){
    ((cout << args<<" "), ...);
    cout << "\n";
}
#define ini(x) int x;cin>>x;
#define inll(x) ll x;cin>>x;
template<typename U>
istream &operator>>(istream &in,vector<U> &v){
    for(auto &i:v) in>>i;
    return in;
}
template<typename U>
ostream &operator<<(ostream &out,vector<U> &v){
    for(auto &i:v) out<<i<<" ";
    return out;
}
#define ins(x) string x;cin>>x;
#define invi(x,n) vi x(n);ff(i,0,n-1) cin>>x[i];
#define invll(x,n) vll x(n);ff(i,0,n-1) cin>>x[i];
// loop snippets
#define ff(i,init,fin) for(int i=init;i<=fin;i++)
#define fb(i,init,fin) for(int i=init;i>=fin;i--)
#define ffs(i,init,fin,step) for(int i=init;i<=fin;i=i+step)
#define fbs(i,init,fin,step) for(int i=init;i>=fin;i=i-step)
#define ffit(it,x) for(auto it=x.begin();it!=x.end();it++)
#define ffa(it,x) for(auto &it:x)
//bit snippets
#define popcnt __builtin_popcount
//function snippets
// ll minimum(ll a,ll b){if(a<b) return a;else return b;}
// ll maximum(ll a,ll b){if(a>b) return a;else return b;}
// ll absolute(ll a){if(a>=0)return a;else return a*-1;}
// ll lcm (ll a, ll b) {return a / gcd(a, b) * b;}
ll mod = 1e9+7;
ll binpow(ll x, ll y,ll M)
{
    if (y == 0)
        return 1;
    ll p = binpow(x, y / 2, M) % M;
    p = (p * p) % M;
    return (y % 2 == 0) ? p : (x * p) % M;
}
ll mod_sub(ll a,ll b){return ((a-b)%mod + mod) % mod;}
ll mod_add(ll a,ll b){return ((a%mod) + (b%mod))%mod;}
ll mod_mult(ll a,ll b){a = a%mod;b=b%mod; return ((a*b)%mod + mod)%mod;}
ll mod_inverse(ll A, ll M)
{
    return binpow(A, M - 2, M);
}
ll mod_divide(ll a, ll b)
{

    a = a % mod;
    ll inv = mod_inverse(b, mod);
    if (inv == -1)
        return 0;
    else
        return mod_mult(inv,a);
}
class mint {
    ll a;

  public:
    mint() { this->a = 0; }
    mint(ll a) { this->a = a; }
    ll &operator()() { return a; }
    mint operator=(ll x) {
        this->a = x;
        return *this;
    }
    mint operator%(ll mod) {
        mint res = *this;
        if (res.a > mod)
            res.a = res.a % mod;
        return res;
    }
    mint operator+(mint b) {
        mint x = (*this) % mod;
        mint y = (b) % mod;
        mint z = (x() + y()) % mod;
        return z = z % mod;
    }
    mint operator-(mint b) {
        return (((*this) % mod)() - (b % mod)() + mod) % mod;
    }
    mint operator*(mint b) {
        mint x = (*this) % mod;
        mint y = (b % mod) % mod;
        mint z = (x() * y());
        return z = z % mod;
    }
    mint operator^(mint b){
        mint x = binpow((*this)(), b(), mod);
        return x;
    }
    mint operator/(mint b){
        mint x = mod_inverse(b(), mod);
        return x = (*this) * x;
    }
    mint operator+=(mint b) { return *this = (*this) + b; }
    mint operator-=(mint b) { return *this = (*this) - b; }
    mint operator*=(mint b) { return *this = (*this) * b; }
    mint operator%=(ll b) { return *this = (*this) % b; }
    mint operator/=(mint b) { return *this = (*this) / b; }
};
typedef vector<mint> vm;
typedef vector<vector<mint>> vvm;
ld prec = 1e-7;
bool iseq(ld v1,ld v2){ return abs(v2 - v1)<= prec; }
bool islt(ld v1,ld v2){if(iseq(v1,v2)) return 0; return v1<v2;}
bool isgt(ld v1,ld v2){if(iseq(v1,v2)) return 0; return v1>v2;}
bool isgte(ld v1,ld v2){if(iseq(v1,v2)) return 1; return v1>v2;}
bool islte(ld v1,ld v2){if(iseq(v1,v2)) return 1; return v1<v2;}
template<typename T,typename U>
U slicing(T const& v,int X, int Y){auto first = v.begin() + X;auto last = v.begin() + Y + 1;auto cont=U(first, last);return cont;}
#define substring slicing<string,string>
#define trace1d(arr,n) cout<<#arr<<"\n";for(int i=0;i<=n;i++)cout<<(arr[i])<<" ";cout<<"\n";
#define trace1d(arr,n) cout<<#arr<<"\n";for(int i=0;i<=n;i++)cout<<(arr[i])<<" ";cout<<"\n";
#define trace2d(arr,n,m) cout<<#arr<<"\n";for(int i=0;i<=n;i++){for(int j=0;j<=m;j++){cout<<(arr[i][j])<<" ";}cout<<"\n";}
#define trace(x) cout<<#x<<" "<<x<<"\n";
#define init1d(val,arr,n) for(int i=0;i<=n;i++){arr[i]=val;}
#define init2d(val,arr,n,m) for(int i=0;i<=n;i++){for(int j=0;j<=m;j++){arr[i][j]=val;}}
template <class T> using oset = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <class T,class U> using omap = tree<T, U, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <class T> using omm = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;
#define ook order_of_key
#define fbo find_by_order
template <typename F>
void deb(F&& lamda){
    #ifndef ONLINE_JUDGE
        lamda();
    #endif
}


/*
 * Questions to ask before submitting any code on OJ
 * Q1. Is my approach handling all the cases ? Think of some edge cases
 * Q2. How complicated is my approach
 * Q3. Will your implementation be a barrier?
 * Remember:
 * Competition is with yourself
 * int overflow, array bounds
 * special cases (n=1?)
 * do smth instead of nothing and stay organized
 * WRITE STUFF DOWN
 * DON'T GET STUCK ON ONE APPROACH
 * No OJ on kickstart
 */
void solve();
int main(){
    #ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);freopen("output.txt", "w", stdout);
    #endif
    fastIO;
    int t;
    int ctr=1;
    cin>>t;
    while(t--){
        deb(
            [&]{
                cout<<"Case #"<<ctr<<" : \n";
                cout.flush();
                ctr++;
            }
        );
        solve();
    }
    deb([]{
        cerr << "Runtime is: " << (clock() * 1.0 / CLOCKS_PER_SEC)*1000 << "ms\n";
    });
}
ll CRT(vpll divisor_remainder){
	ll n = divisor_remainder.size();
	ll m = 1;
	ff(i,0,n-1){
		m*=divisor_remainder[i].fi;
	}
	ll ans=0;
	ff(i,0,n-1){
		ll curr = 1;
		ff(j,0,i-1){
            curr*=divisor_remainder[j].fi;
			curr *=mod_inverse(divisor_remainder[j].first,divisor_remainder[i].first);
			curr%=m;
		}
		ans+= ( (divisor_remainder[i].se - ans + m)*curr)%m;
		ans%=m;
	}
	return ans;
}
class DSU {
  public:
    vi Rank, Par;
    vvi Ele;
    DSU(int n) {
        Par.resize(n + 1);
        Rank.resize(n + 1);
        for (int i = 1; i <= n; i++) {
            Par[i] = i;
            Rank[i] = 1;
        }
    }
    int find_leader(int);
    void union_sets(int, int);
};
int DSU::find_leader(int v) {
    if (Par[v] == v) {
        return v;
    }
    return Par[v] = find_leader(Par[v]);
}
void DSU::union_sets(int a, int b) {
    a = find_leader(a);
    b = find_leader(b);
    if (a == b)
        return;
    if (Rank[a] < Rank[b])
        swap(a, b);
    Par[b] = a;
    if (Rank[a] == Rank[b])
        Rank[a]++;
}

template<typename T>
struct Matrix{
    vector<vector<T>> a;
    // multiplication of matrices
    Matrix<T> operator *(Matrix<T> &vec){
        auto b = vec.a;
        int a_r = a.size(),a_c = a[0].size();
        int b_r = b.size(),b_c = b[0].size();

        Matrix<T> res{vector<vector<T>>(a_r,vector<T>(b_c))};
        ff(i,0,a_r-1){
            ff(j,0,b_c-1){
                ff(k,0,b_c-1){
                    res.a[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return res;
    }

    // addition
    Matrix<T> operator +(Matrix<T> &vec){
        auto b = vec.a;
        int a_r = a.size(),a_c = a[0].size();
        Matrix<T> res{vector<vector<T>>(a_r,vector<T>(a_c))};
        ff(i,0,a_r-1){
            ff(j,0,a_c-1){
                res.a[i][j] = a[i][j] + b[i][j];
            }
        }
        return res;
    }

    // substraction
    Matrix<T> operator -(Matrix<T> &vec){
        auto b = vec.a;
        int a_r = a.size(),a_c = a[0].size();
        Matrix<T> res{vector<vector<T>>(a_r,vector<T>(a_c))};
        ff(i,0,a_r-1){
            ff(j,0,a_c-1){
                res.a[i][j] = a[i][j] - b[i][j];
            }
        }
        return res;
    }

};
#define root idx
#define lc 2*idx+1
#define rc 2*idx+2
class segtree {
  public:
    vll stree, arr;
    int n;

    void update(int, int, int, int, int);
    void build(int, int, int);
    ll query(int, int, int, int, int);

    segtree() {}
    segtree(vll _arr) {
        n = _arr.size();
        arr = _arr;
        stree.resize(4 * n);
        build(0, n - 1, 0);
    }
};

void segtree::update(int i, int val, int l, int r, int idx) {
    if (l == r) {
        if (l == i) {
            stree[root] = arr[i] = val;
        }
        return;
    }
    if (i > r || i < l)
        return;
    int mid = (l + r) / 2;
    update(i, val, l, mid, lc);
    update(i, val, mid + 1, r, rc);

    // implement how two segments are begin merged
    // stree[root] =
}

void segtree::build(int l, int r, int idx) {
    if (l == r) {
        stree[root] = arr[l];
    }
    if (l < r) {
        int mid = (l + r) / 2;
        build(l, mid, lc);
        build(mid + 1, r, rc);

        // implement how two segments are begin merged
    }
}

ll segtree::query(int l, int r, int idx, int curr_l, int curr_r) {
    if (curr_l > r || curr_r < l) {
        return 0;
    } else if (l <= curr_l && curr_r <= r) {
        return stree[root];
    } else {
        int mid = (curr_l + curr_r) / 2;

        // implement merging
    }
}
vector<bool> make_sieve(ll n){
	//sieve version 1 TC --> nloglogn the impact of the constant factor is high in this version
	vector<bool> is_prime(n+1, true);
	is_prime[0] = is_prime[1] = false;
	//for (int i = 2; i <= n; i++) {
		//if (is_prime[i] && (long long)i * i <= n) {
			//for (int j = i * i; j <= n; j += i)
				//is_prime[j] = false;
		//}
	//}
	for (int i = 2; i * i <= n; i++) {
		if (is_prime[i]) {
	        for (int j = i * i; j <= n; j += i)
			    is_prime[j] = false;
		}
	}
	return is_prime;
}

vector<bool> segmentedSieve(long long L, long long R) {
    // generate all primes up to sqrt(R)
    long long lim = sqrt(R);
    vector<bool> mark(lim + 1, false);
    vector<long long> primes;
    for (long long i = 2; i <= lim; ++i) {
        if (!mark[i]) {
            primes.emplace_back(i);
            for (long long j = i * i; j <= lim; j += i)
                mark[j] = true;
        }
    }

    vector<bool> isPrime(R - L + 1, true);
    for (long long i : primes)
        for (long long j = max(i * i, (L + i - 1) / i * i); j <= R; j += i)
            isPrime[j - L] = false;
    if (L == 1)
        isPrime[0] = false;
    return isPrime;
}
ld SQRT(ll x){
	if(x==0) return 0;
	ld l = prec,r = x/2+1,res=0;
	while(isgt(r-l,prec)){
		ld mid = (r+l)/2;
		if(islte(mid*mid,x)){
			l=mid;
			res = mid;
		}
		else{
			r=mid;
		}
	}
	return res;
}
ll SQRTL(ll x){
	if(x==0) return 0;
    	ll l = 1,r = x/2+1;
    	while(r-l>1){
    		ld mid = (r+l)/2;
    		if(mid*mid<=x){
    			l=mid;
    		}
    		else{
    			r=mid;
    		}
    	}
    	if(r-l<=1){
	        if(r*r<=x){
    	        return r;
            }
	    }
        return l;
}
void solve(){

}