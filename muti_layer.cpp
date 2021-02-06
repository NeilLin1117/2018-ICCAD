#include <iostream>
#include <list>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <math.h>
#include <sstream>
#include <fstream>
#include <limits.h>
#include<map>
#include <time.h>
using namespace std;
int global_count = 0;


struct Edge {
	int start, end;
	int start_serial, end_serial;
	double length;
	bool visit;
};

vector<Edge> edge;
struct Points {
	double x, y;
};



vector<int> father; int  Count;

class grid
{
public:
	grid() :value(0), visited(false), conbine(false) {}
	int i;
	int j;
	int k;
	bool visited;
	int value;
	grid(int i, int j, int k) {
		value = 0;
		visited = false;
		this->i = i;
		this->j = j;
		this->k = k;
		//this->conbine = conbine;

	}
	int fi;
	int fj;
	int fk;
	double x_point;
	double y_point;
	int sym;
	bool conbine;
};


class pin {
public:
	pin() :x(0), y(0), serial_number(0), metal_layer(0) {}
	pin &operator = (const pin &inpin) {
		x = inpin.x;
		y = inpin.y;
		serial_number = inpin.serial_number;
		metal_layer = inpin.metal_layer;
	}
	int serial_number;
	int metal_layer;
	double x;
	double y;
};

class net {
public:
	net() :priority(0) {}
	int serial_number;
	vector <pin> connect;
	char critical;
	int priority;
	bool operator<(const net &rhs) const {
		if ((critical == 'Y' && rhs.critical == 'Y') ||
			(critical == 'N' && rhs.critical == 'N'))
			return connect.size() > rhs.connect.size();
		else {
			if (critical == 'Y')
				return true;
			else
				return false;
		}
	}
};

class mul_net {
public:
	mul_net() :priority(0) {}
	int serial_number;
	vector <pin> connect;
	char critical;
	int priority;
	vector<Edge> edges;
	bool operator<(const mul_net &rhs) const {
		if ((critical == 'Y' && rhs.critical == 'Y') ||
			(critical == 'N' && rhs.critical == 'N'))
			return connect.size() > rhs.connect.size();
		else {
			if (critical == 'Y')
				return true;
			else
				return false;
		}
	}
};


class blockage {
public:
	blockage() :x1(0), y1(0), x2(0), y2(0) {}
	double x1;
	double y1;
	double x2;
	double y2;
};

vector <string> string_data;
vector <pin> pin_Data;
vector <net> net_Data;
vector <mul_net> multi_pin_net;
vector <blockage> blk_Data;
//ofstream out;
int max_x = 0;
int max_y = 0;

grid ***layer_ptr = new grid **[4];
vector <grid> find_neighbor(grid point);
void traceback(grid s, grid t);
int count_words(string s);
int count_field_num(pin p1, pin p2);
void init_value(int net, char pinset,int serial);
void print();
vector<vector <int> > rand_order(int n, int times);


////////////////////////////////////Cost Function//////////////////////////////////
struct clear_net {
	Points start;
	Points end;
	int layer;
	int color;
};

map<int, vector<clear_net> > net_of_storage;
map<int, vector<clear_net> > copy_net_of_storage;
map<int, vector<clear_net> >::iterator iter;

double Detour_critical = 0.0;
int global_flag = 1;
const int MAXN = 5005;
double total_wire_length = 0.0;
double total_via = 0.0;
double balance_score = 0.0;
bool find_error = false;


class balance_wire {
public:
	balance_wire() :total(0.0), color1(0.0), color2(0.0) {};
	double total;
	double color1;
	double color2;
	double ratio1;
	double ratio2;
};

struct Cost_Edge {            
	int start, end;      
	double length;  
	double wire;      
	double via;
	bool visit;         
};

vector<Cost_Edge> edge_data;
map<int, balance_wire > balance_ratio;
map<int, balance_wire > ::iterator iter_cost;

struct Cost_Point {           
	double x, y, z;        
	int color;
	int flag;               
};
vector<Cost_Point> points_data;	

int cost_father[MAXN];      
							

double  getPriority(Cost_Point a, Cost_Point b)  
{
	double len;
	//len = sqrt(double((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)));
	len = double(abs(a.x - b.x) + abs(a.y - b.y) + abs(a.z - b.z));
	return len;
}

double  get_wire_length(Cost_Point a, Cost_Point b)  
{
	double len;
	len = double(abs(a.x - b.x) + abs(a.y - b.y));
	return len;
}

double  getVia(Cost_Point a, Cost_Point b)  
{
	double len;
	len = double(abs(a.z - b.z));
	return len;
}

bool cost_cmp(Cost_Edge a, Cost_Edge b) {      
	return a.length<b.length;
}

int cost_find(int x) {             
	if (cost_father[x] == x) return x;
	cost_father[x] = cost_find(cost_father[x]);
	return cost_father[x];
}

bool cost_Union(int x, int y)
{
	int f1 = cost_find(x);  
	int f2 = cost_find(y);  
	if (f1 == f2) return false;  
	else if (f1<f2) cost_father[f1] = f2;  
	else cost_father[f2] = f1;
	return true;
}

void cost_kruskal(int n) { 
	int i, j = 0;
	double sum = 0;
	for (i = 0; i<n; i++) cost_father[i] = i;   
	sort(edge_data.begin(), edge_data.end(), cost_cmp);      
	for (i = 0; i<edge_data.size() && j<n; i++)            
	{
		if (cost_Union(edge_data[i].start, edge_data[i].end)) {

			edge_data[i].visit = 1;      
			j++;
		}
	}
	//return sum;                    
}

void setting(int start, int end) {
	points_data[start].flag = 1;
	points_data[end].flag = 1;
	Cost_Edge edge;
	edge.start = start;
	edge.end = end;
	edge.length = getPriority(points_data[start], points_data[end]);
	edge.wire = get_wire_length(points_data[start], points_data[end]);
	edge.via = getVia(points_data[start], points_data[end]);
	edge.visit = 0;
	edge_data.push_back(edge);
}

void setEdge(clear_net data) {
	if (data.start.x == data.end.x && data.start.y == data.end.y) {
		Cost_Point point;
		point.x = data.start.x;
		point.y = data.start.y;
		point.z = data.layer;
		point.color = data.color;
		point.flag = 0;
		points_data.push_back(point);
		//P++;
	}
	else {
		Cost_Point point;
		point.x = data.start.x;
		point.y = data.start.y;
		point.z = data.layer;
		point.color = data.color;
		point.flag = 1;
		points_data.push_back(point);
		//P++;
		point.x = data.end.x;
		point.y = data.end.y;
		point.z = data.layer;
		point.color = data.color;
		point.flag = 1;
		points_data.push_back(point);
		setting(points_data.size() - 2, points_data.size() - 1);
	}
}

void set_steiner(int start, int ch_start, int ch_end) {
	if (points_data[start].x == points_data[ch_start].x &&
		((points_data[start].y <  points_data[ch_start].y &&  points_data[start].y >  points_data[ch_end].y)
			|| (points_data[start].y >  points_data[ch_start].y &&  points_data[start].y <  points_data[ch_end].y))
		&& abs(points_data[start].z - points_data[ch_start].z) == 1) {
		Cost_Point point;
		point.x = points_data[start].x;
		point.y = points_data[start].y;
		point.z = points_data[ch_start].z;
		point.color = points_data[ch_start].color;
		points_data.push_back(point);
		setting(start, points_data.size() - 1);
		setting(points_data.size() - 1, ch_start);
		setting(points_data.size() - 1, ch_end);

	}

	if (points_data[start].y == points_data[ch_start].y &&
		((points_data[start].x <  points_data[ch_start].x &&  points_data[start].x >  points_data[ch_end].x)
			|| (points_data[start].x >  points_data[ch_start].x &&  points_data[start].x <  points_data[ch_end].x))
		&& abs(points_data[start].z - points_data[ch_start].z) == 1) {
		Cost_Point point;
		point.x = points_data[start].x;
		point.y = points_data[start].y;
		point.z = points_data[ch_start].z;
		point.color = points_data[ch_start].color;
		points_data.push_back(point);
		setting(start, points_data.size() - 1);
		setting(points_data.size() - 1, ch_start);
		setting(points_data.size() - 1, ch_end);

	}
}

void count_balance() {
	for (int i = 0; i < edge_data.size(); i++) {
		if (edge_data[i].visit) {
			int start = edge_data[i].start;
			int end = edge_data[i].end;
			if (points_data[start].z == points_data[end].z &&  points_data[start].color == points_data[end].color) {
				int layer = (int)points_data[start].z;
				balance_ratio[layer].total += edge_data[i].wire;
				if (points_data[start].color == 1) {
					balance_ratio[layer].color1 += edge_data[i].wire;
				}
				else {
					balance_ratio[layer].color2 += edge_data[i].wire;
				}
			}
		}
	}
}

void cost_pre_kruskal(int critical) {
	for (int i = 0; i< points_data.size() - 1; i++)      
		for (int j = i + 1; j < points_data.size(); j++) {
			if (points_data[i].x == points_data[j].x &&
				points_data[i].y == points_data[j].y &&
				(abs(points_data[i].z - points_data[j].z) == 1)) {
				setting(i, j);
			}
		}
	vector<Cost_Edge> copy_data;
	for (int i = 0; i < edge_data.size(); i++) {
		copy_data.push_back(edge_data[i]);
	}

	
	for (int i = 0; i<copy_data.size() - 1; i++) {
		int start = copy_data[i].start;
		int end = copy_data[i].end;
		for (int j = i + 1; j < copy_data.size(); j++) {

			int ch_start = copy_data[j].start;
			int ch_end = copy_data[j].end;
			if ((points_data[start].x != points_data[ch_start].x && points_data[start].y != points_data[ch_start].y
				&& abs(points_data[start].z - points_data[ch_start].z) == 1) && (points_data[start].x != points_data[ch_end].x && points_data[start].y != points_data[ch_end].y
					&& abs(points_data[start].z - points_data[ch_end].z) == 1)) {
			
				if ((points_data[end].x != points_data[ch_start].x && points_data[end].y != points_data[ch_start].y
					&& abs(points_data[end].z - points_data[ch_start].z) == 1) && (points_data[end].x != points_data[ch_end].x
						&& points_data[end].y != points_data[ch_end].y && abs(points_data[end].z - points_data[ch_end].z) == 1))
				{
					
					if (points_data[start].x == points_data[end].x && points_data[ch_start].y == points_data[ch_end].y) {
						if ((points_data[start].x < points_data[ch_start].x && points_data[start].x > points_data[ch_end].x)
							|| (points_data[start].x > points_data[ch_start].x && points_data[start].x < points_data[ch_end].x))
						{
							if ((points_data[ch_start].y < points_data[start].y && points_data[ch_start].y > points_data[end].y)
								|| (points_data[ch_start].y > points_data[start].y && points_data[ch_start].y < points_data[end].y)) {
								Cost_Point point;
								point.x = points_data[start].x;
								point.y = points_data[ch_start].y;
								point.z = points_data[start].z;
								point.color = points_data[start].color;
								points_data.push_back(point);
								setting(start, points_data.size() - 1);
								setting(end, points_data.size() - 1);
								point.x = points_data[start].x;
								point.y = points_data[ch_start].y;
								point.z = points_data[ch_start].z;
								point.color = points_data[ch_start].color;
								points_data.push_back(point);
								setting(points_data.size() - 2, points_data.size() - 1);
								setting(ch_start, points_data.size() - 1);
								setting(ch_end, points_data.size() - 1);
							}
						}
					}
					if (points_data[start].y == points_data[end].y && points_data[ch_start].x == points_data[ch_end].x) {
						if ((points_data[start].y < points_data[ch_start].y && points_data[start].y > points_data[ch_end].y)
							|| (points_data[start].y > points_data[ch_start].y && points_data[start].y < points_data[ch_end].y))
						{
							if ((points_data[ch_start].x < points_data[start].x && points_data[ch_start].x > points_data[end].x)
								|| (points_data[ch_start].x > points_data[start].x && points_data[ch_start].x < points_data[end].x)) {
								Cost_Point point;
								point.x = points_data[ch_start].x;
								point.y = points_data[start].y;
								point.z = points_data[start].z;
								point.color = points_data[start].color;
								points_data.push_back(point);
								setting(start, points_data.size() - 1);
								setting(end, points_data.size() - 1);
								point.x = points_data[ch_start].x;
								point.y = points_data[start].y;
								point.z = points_data[ch_start].z;
								point.color = points_data[ch_start].color;
								points_data.push_back(point);
								setting(points_data.size() - 2, points_data.size() - 1);
								setting(ch_start, points_data.size() - 1);
								setting(ch_end, points_data.size() - 1);
							}
						}
					}
				}
			}
		}
	}

	
	for (int i = 0; i < copy_data.size(); i++) {
		int start = copy_data[i].start;
		int end = copy_data[i].end;
		for (int j = 0; j < copy_data.size(); j++) {
			if (i == j)
				continue;
			int ch_start = copy_data[j].start;
			int ch_end = copy_data[j].end;
			set_steiner(start, ch_start, ch_end);
			set_steiner(end, ch_start, ch_end);
		}
	}
	for (int i = 0; i < points_data.size(); i++) {
		if (points_data[i].flag == 0) {
			for (int j = 0; j < copy_data.size(); j++) {
				int ch_start = copy_data[j].start;
				int ch_end = copy_data[j].end;
				set_steiner(i, ch_start, ch_end);
			}
		}
	}



	//cout << "done1!!!" << endl;
	cost_kruskal(points_data.size());
	//cout << "done2!!!" << endl;
	double Via = 0, wire_length = 0, total = 0;
	for (int i = 0; i < edge_data.size(); i++) {
		if (edge_data[i].visit) {
			Via += edge_data[i].via;
			wire_length += edge_data[i].wire;
			total += edge_data[i].length;
		}
		cost_Union(edge_data[i].start, edge_data[i].end);
	}
	int flag = 1;
	for (int i = 0; i < points_data.size(); i++) {
		if (cost_father[i] != points_data.size() - 1) {
			flag = 0;
			break;
		}
	}
	if (flag) {
	
		total_wire_length += wire_length;
		total_via += Via;
		count_balance();
		if (critical) {
			Detour_critical += wire_length;
		}

	}
	else {
		global_flag = 0;
		
	}
}

void count_balance_ratio() {
	for (iter_cost = balance_ratio.begin(); iter_cost != balance_ratio.end(); iter_cost++) {
		double ratios;
		//cout << "Layer " << iter_cost->first << ":" << endl;
		//cout << "Total: " << iter_cost->second.total << endl;
		//cout << "color1: " << iter_cost->second.color1 << endl;
		//cout << "color2: " << iter_cost->second.color2 << endl;
		iter_cost->second.ratio1 = (double)iter_cost->second.color1 / (double)iter_cost->second.total;
		iter_cost->second.ratio2 = (double)iter_cost->second.color2 / (double)iter_cost->second.total;
		//cout << "color1 ratio: " << fixed << setprecision(2) << iter_cost->second.ratio1 << endl;
		//cout << "color2 ratio: " << fixed << setprecision(2) << iter_cost->second.ratio2 << endl << endl;
		ratios = 100 - (abs(iter_cost->second.ratio1 - iter_cost->second.ratio2) * 100);
		balance_score += ratios;
	}
}

void vertify() {
	for (int i = 0; i < multi_pin_net.size(); i++) {
		int serial = multi_pin_net[i].serial_number;
		//cout << "Net: " << serial << endl;
		edge_data.clear();
		points_data.clear();
		for (int j = 0; j < net_of_storage[serial].size(); j++) {
			setEdge(net_of_storage[serial][j]);
		}
		if (multi_pin_net[i].critical == 'Y')
			cost_pre_kruskal(1);
		else
			cost_pre_kruskal(0);
	}

	//cout << "sigal" << endl;

	for (int i = 0; i < net_Data.size(); i++) {
		//cout << "in_net: " << net_Data[i].serial_number << endl;
		int serial = net_Data[i].serial_number;
		//cout << "Net: " << serial << endl;
		edge_data.clear();
		points_data.clear();


		for (int j = 0; j < net_of_storage[serial].size(); j++) {

			setEdge(net_of_storage[serial][j]);
		}
		if (net_Data[i].critical == 'Y')
			cost_pre_kruskal(1);
		else
			cost_pre_kruskal(0);
	}

}


////////////////////////////////////Cost Function//////////////////////////////////

bool cmp(Edge a, Edge b) {
	return a.length<b.length;
}

double getlength(pin a, pin b)
{
	double len;
	len = sqrt(double((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)));
	return len;
}

int find(int x) {
	if (father[x] == x) return x;
	father[x] = find(father[x]);
	return father[x];
}

bool Union(int x, int y)
{
	int f1 = find(x);
	int f2 = find(y);
	if (f1 == f2) return false;
	else if (f1<f2) father[f1] = f2;
	else father[f2] = f1;
	return true;
}

void kruskal(int n) {
	int i, j = 0;
	double sum = 0;
	for (i = 0; i<n; i++) father[i] = i;
	sort(edge.begin(), edge.end(), cmp);
	for (i = 0; i<Count &&j<n; i++)
	{
		if (Union(edge[i].start, edge[i].end)) {

			edge[i].visit = 1;
			j++;
		}
	}
}

void init_kruskal(mul_net &a) {
	int P = a.connect.size();
	int count = 0;
	edge.clear();
	father.clear();
	father.resize(P);
	for (int i = 0; i<P - 1; i++)
		for (int j = i + 1; j < P; j++) {
			count++;
		}
	edge.resize(count);
}

void pre_kruskal(mul_net &a) {
	a.edges.clear();
	init_kruskal(a);
	Count = 0;
	int P = a.connect.size();
	for (int i = 0; i<P - 1; i++)
		for (int j = i + 1; j < P; j++)
		{
			edge[Count].start = i;  edge[Count].end = j;
			edge[Count].start_serial = a.connect[i].serial_number;
			edge[Count].end_serial = a.connect[j].serial_number;
			edge[Count].length = getlength(a.connect[i], a.connect[j]);
			edge[Count].visit = 0;
			Count++;
		}
	kruskal(P);
	for (int i = 0; i < Count; i++) {
		if (edge[i].visit) {
			a.edges.push_back(edge[i]);
		}
	}
}

int Pioneer(grid src)
{
	list<grid>queue;
	layer_ptr[src.i][src.j][src.k].visited = true;
	queue.push_back(src);

	int or_i = src.i;
	int or_j = src.j;
	int or_k = src.k;

	grid sta = src;

	int tmp_sym = layer_ptr[src.i][src.j][src.k].sym;

	while (!queue.empty()) {

		src = queue.front();
		queue.pop_front();
		if (

			//(src.i == tgt.i && src.j == tgt.j && src.k == tgt.k)

			layer_ptr[src.i][src.j][src.k].sym == tmp_sym &&
			!(src.i == sta.i && src.j == sta.j && src.k == sta.k) &&
			layer_ptr[src.i][src.j][src.k].conbine == true
			) {
			//cout << "find point " << layer_ptr[src.i][src.j][src.k].sym << endl;
			global_count++;
			src.sym = tmp_sym;
			//cout << "src.sym  " << src.sym << endl;
			traceback(sta, src);
			return layer_ptr[src.i][src.j][src.k].value;
			break;
		}

		vector <grid> neighbor = find_neighbor(src);
		for (int i = 0; i<neighbor.size(); i++) {
			int li = neighbor[i].i;
			int lj = neighbor[i].j;
			int lk = neighbor[i].k;

			if (layer_ptr[li][lj][lk].visited == false &&
				(
					layer_ptr[li][lj][lk].sym == -2 ||
					layer_ptr[li][lj][lk].sym == -1 ||
					layer_ptr[li][lj][lk].sym == tmp_sym
					)
				) {
				layer_ptr[li][lj][lk].visited = true;
				layer_ptr[li][lj][lk].value = src.value + 1;
				layer_ptr[li][lj][lk].fi = src.i;
				layer_ptr[li][lj][lk].fj = src.j;
				layer_ptr[li][lj][lk].fk = src.k;
				queue.push_back(layer_ptr[li][lj][lk]);

			}
		}
	}
	//////////////////////////////////////////////// reroute area ///////////////////////////////////////////////
	find_error = true;
	//getchar();
}

void traceback(grid s, grid t)
{
	int tmpi = t.i;
	int tmpj = t.j;
	int tmpk = t.k;
	int net = layer_ptr[t.i][t.j][t.k].sym;
	//cout << "in_net _traceback " << endl;
	while (tmpi != s.i || tmpj != s.j || tmpk != s.k) {
		layer_ptr[t.i][t.j][t.k].sym = net;
		layer_ptr[t.i][t.j][t.k].conbine = true;
		tmpi = layer_ptr[t.i][t.j][t.k].fi;
		tmpj = layer_ptr[t.i][t.j][t.k].fj;
		tmpk = layer_ptr[t.i][t.j][t.k].fk;
		t.i = tmpi;
		t.j = tmpj;
		t.k = tmpk;
	}

	layer_ptr[s.i][s.j][s.k].sym = net;
	layer_ptr[t.i][t.j][t.k].conbine = true;
}

vector <grid> find_neighbor(grid point)
{
	vector <grid> neighbor;

	if (point.i % 2 == 0) {

		if (point.k + 1 < max_x + 1) {
			if (layer_ptr[point.i][point.j][point.k + 1].sym != -3 && layer_ptr[point.i][point.j][point.k + 1].sym != -4)
				neighbor.push_back(layer_ptr[point.i][point.j][point.k + 1]);
		}
		if (point.k - 1 >= 0) {

			if (layer_ptr[point.i][point.j][point.k - 1].sym != -3 &&layer_ptr[point.i][point.j][point.k - 1].sym != -4)
				neighbor.push_back(layer_ptr[point.i][point.j][point.k - 1]);
		}

	}
	else {

		if (point.j + 1 < max_y + 1) {

			if (layer_ptr[point.i][point.j + 1][point.k].sym != -3&&layer_ptr[point.i][point.j + 1][point.k].sym != -4)
				neighbor.push_back(layer_ptr[point.i][point.j + 1][point.k]);
		}
		if (point.j - 1 >= 0) {

			if (layer_ptr[point.i][point.j - 1][point.k].sym != -3&&layer_ptr[point.i][point.j - 1][point.k].sym != -4)
				neighbor.push_back(layer_ptr[point.i][point.j - 1][point.k]);
		}

	}
	if (point.i + 1 < 4) {
		if (layer_ptr[point.i + 1][point.j][point.k].sym != -3 &&layer_ptr[point.i + 1][point.j][point.k].sym != -4) {

			neighbor.push_back(layer_ptr[point.i + 1][point.j][point.k]);
		}
	}
	if (point.i - 1 >= 0) {
		if (layer_ptr[point.i - 1][point.j][point.k].sym != -3 &&layer_ptr[point.i - 1][point.j][point.k].sym != -4) {

			neighbor.push_back(layer_ptr[point.i - 1][point.j][point.k]);
		}
	}

	if (neighbor.empty()) {
		cout << "empty" << endl;
		find_error = true;
	}

	return neighbor;
}

void Parser(int choice)
{
	if (choice == 0 || choice == 2) {

		for (int i = 0; i < string_data.size(); i++) {
			stringstream ss(string_data[i]);
			if (choice == 0) {
				pin a;
				ss >> a.serial_number;
				ss >> a.metal_layer;
				ss >> a.x;
				ss >> a.y;
				if (a.x > max_x)
					max_x = a.x;
				if (a.y > max_y)
					max_y = a.y;
				pin_Data.push_back(a);
			}
			else {
				//stringstream ss(string_data[i]);
				blockage a;
				ss >> a.x1;
				ss >> a.y1;
				ss >> a.x2;
				ss >> a.y2;
				blk_Data.push_back(a);
			}
		}
	}
	else {
		for (int i = 0; i < string_data.size(); i++) {
			int str_size = count_words(string_data[i]);
			//cout << str_size << endl;
			stringstream ss(string_data[i]);
			if (str_size == 3) {
				net a;
				ss >> a.serial_number;
				int tmp;
				for (int i = 0; i < str_size - 1; i++) {
					ss >> tmp;
					a.connect.push_back(pin_Data[tmp - 1]);
				}
				ss >> a.critical;
				for (int j = 0; j < a.connect.size() - 1; j++) {

					int sub_x = abs(a.connect[j].x - a.connect[j + 1].x);
					int sub_y = abs(a.connect[j].y - a.connect[j + 1].y);
					a.priority += sub_x + sub_y;
				}
				net_Data.push_back(a);
			}
			else
			{
				mul_net a;
				ss >> a.serial_number;
				int tmp;
				for (int i = 0; i < str_size - 1; i++) {
					ss >> tmp;
					a.connect.push_back(pin_Data[tmp - 1]);
				}
				ss >> a.critical;
				for (int j = 0; j < a.connect.size() - 1; j++) {

					int sub_x = abs(a.connect[j].x - a.connect[j + 1].x);
					int sub_y = abs(a.connect[j].y - a.connect[j + 1].y);
					a.priority += sub_x + sub_y;
				}
				multi_pin_net.push_back(a);
			}
		}
	}
}


int count_words(string s)
{
	int words = 0;
	for (int i = 0; i < s.length(); i++) {
		if (s[i] == ' ') {
			words++;
		}
	}
	return words;
}

void readFile(const char* file, int choice)
{
	string input;
	ifstream infile(file, ios::in);
	if (!infile) // overloaded ! operator
	{
		cerr << "File could not be opened" << endl;
		exit(EXIT_FAILURE);
	}
	string_data.clear();
	string tmp;
	while (getline(infile, tmp)) {
		if (infile.eof()) {
			if (tmp != "")
				string_data.push_back(tmp);
			break;
		}

		string_data.push_back(tmp);
	}
	infile.close();
	Parser(choice);
}


void find_grid(int sym,  int Serial) {

	//cout << "in_find_grad" << Serial << endl;
	vector<clear_net> tmp;
	for (int i = 0; i < 4; i++) {
		if (i % 2 == 0) {
			for (int j = max_y; j >= 0; j--) {

				int color = j % 2;
				int start = 0, end = 0;
				for (; start < max_x + 1; start++) {
					clear_net copy_tmp;
					if (layer_ptr[i][j][start].sym == sym) {

						end = start;
						copy_tmp.start.x = start * 0.5;
						copy_tmp.start.y = j * 0.5;

						while (layer_ptr[i][j][end].sym == sym && end < max_x + 1)
							end++;
						copy_tmp.end.x = (end - 1)*0.5;
						copy_tmp.end.y = j * 0.5;
						copy_tmp.layer = i + 1;
						copy_tmp.color = color + 1;

						if (start - 1 >= 0)
							layer_ptr[i][j][start - 1].sym = -4;

						if (end  < max_x + 1)
							layer_ptr[i][j][end].sym = -4;
						start = end;
						tmp.push_back(copy_tmp);
					}

				}
			}
		}
		else {
			for (int k = 0; k < max_x + 1; k++) {
				int start = max_y, end = 0;
				int color = k % 2;
				for (; start >= 0; start--) {
					clear_net copy_tmp;
					if (layer_ptr[i][start][k].sym == sym) {
						end = start;
						copy_tmp.start.x = k * 0.5;
						copy_tmp.start.y = start * 0.5;

						while (layer_ptr[i][end][k].sym == sym && end >= 0)
							end--;
						copy_tmp.end.x = k * 0.5;
						copy_tmp.end.y = (end + 1) * 0.5;
						copy_tmp.layer = i + 1;
						copy_tmp.color = color + 1;

						if (end >= 0)
							layer_ptr[i][end][k].sym = -4;
						if (start + 1 < max_y + 1)
							layer_ptr[i][start + 1][k].sym = -4;
						start = end;
						tmp.push_back(copy_tmp);
					}

				}
			}
		}
	}
	net_of_storage.insert(pair<int, vector<clear_net> >(Serial, tmp));

}


void print_net_of_storage(ofstream & out ,int score,char* file) {
	out << file << endl;
	for (iter = copy_net_of_storage.begin(); iter != copy_net_of_storage.end(); iter++) {
		out << "Net " << iter->first << endl;
		for (int j = 0; j<iter->second.size(); j++) {


			out << fixed << setprecision(1) << iter->second[j].start.x << ' ' << iter->second[j].start.y << ' '
				<< iter->second[j].end.x << ' ' << iter->second[j].end.y << ' '
				<< iter->second[j].layer << ' ' << iter->second[j].color << endl;

		}
	}
	out.close();
}

void print() {
	for (int i = 0; i<4; i++) {
		for (int j = max_y; j >= 0; j--) {
			for (int k = 0; k<max_x + 1; k++) {
				cout << layer_ptr[i][j][k].sym;
			}
			cout << endl;
		}
		cout << 'l' << endl;
	}
}

void init_circuit()
{
	for (int i = 0; i < 4; i++) {
		layer_ptr[i] = new grid*[max_y + 1];
		for (int j = 0; j < max_y + 1; j++)
			layer_ptr[i][j] = new grid[max_x + 1]();
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < max_y + 1; j++) {
			for (int k = 0; k < max_x + 1; k++) {
				layer_ptr[i][j][k].x_point = k * 0.5;
				layer_ptr[i][j][k].y_point = j * 0.5;
			}
		}
	}

	for (int i = 0; i<4; i++) {
		if (i % 2 == 0) {
			for (int j = 0; j<max_y + 1; j++) {
				for (int k = 0; k<max_x + 1; k++) {
					layer_ptr[i][j][k].sym = -1;
					layer_ptr[i][j][k].i = i;
					layer_ptr[i][j][k].j = j;
					layer_ptr[i][j][k].k = k;
				}
			}
		}
		else {
			for (int j = 0; j<max_y + 1; j++) {
				for (int k = 0; k<max_x + 1; k++) {
					layer_ptr[i][j][k].sym = -2;
					layer_ptr[i][j][k].i = i;
					layer_ptr[i][j][k].j = j;
					layer_ptr[i][j][k].k = k;
				}
			}
		}
	}


	for (int t = 0; t < blk_Data.size(); t++) {
		for (int j = (blk_Data[t].y1 * 2); j <= (blk_Data[t].y2 * 2); j++) {
			for (int k = (blk_Data[t].x1 * 2); k <= (blk_Data[t].x2 * 2); k++)
				for (int i = 0; i<4; i++)
					layer_ptr[i][j][k].sym = -3;
		}
	}

	for (int i = 0; i<net_Data.size(); i++) {
		for (int j = 0; j<net_Data[i].connect.size(); j++) {
			pin tmp = net_Data[i].connect[j];
			int y = tmp.y * 2;
			int x = tmp.x * 2;
			layer_ptr[tmp.metal_layer - 1][y][x].sym = net_Data[i].serial_number;

		}
	}

	for (int i = 0; i<multi_pin_net.size(); i++) {
		for (int j = 0; j<multi_pin_net[i].connect.size(); j++) {
			pin tmp = multi_pin_net[i].connect[j];
			int y = tmp.y * 2;
			int x = tmp.x * 2;
			layer_ptr[tmp.metal_layer - 1][y][x].sym = multi_pin_net[i].serial_number;
		}
	}

	sort(net_Data.begin(), net_Data.end());
	sort(multi_pin_net.begin(), multi_pin_net.end());
	for (int i = 0; i < multi_pin_net.size(); i++) {

		pre_kruskal(multi_pin_net[i]);

	}

}

void init_travelsal()
{
	for (int i = 0; i<4; i++) {
		for (int j = max_y; j >= 0; j--) {
			for (int k = 0; k<max_x + 1; k++) {
				layer_ptr[i][j][k].value = 0;
				layer_ptr[i][j][k].visited = false;
				layer_ptr[i][j][k].fi = 0;
				layer_ptr[i][j][k].fj = 0;
				layer_ptr[i][j][k].fk = 0;
			}
			//cout << endl;
		}
		//cout << 'l' << endl;
	}
}


void init_value(int net, char pinset, int serial)
{
	int serials = serial;
	init_travelsal();
	if (pinset == 'M') {
		for (int j = 0; j < multi_pin_net[serials].connect.size(); j++) {
			pin tmp = multi_pin_net[serials].connect[j];
			int y = tmp.y * 2;
			int x = tmp.x * 2;
			if (tmp.metal_layer - 1 == 0 || tmp.metal_layer - 1 == 2) {

				if (x + 1 < max_x + 1)
				{
					if (layer_ptr[tmp.metal_layer - 1][y][x + 1].sym == -3) {

						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -1;
					}

					else if (layer_ptr[tmp.metal_layer - 1][y][x + 1].sym == -4) {

						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -3;
					}
				}
				if (x - 1 >= 0) {
					if (layer_ptr[tmp.metal_layer - 1][y][x - 1].sym == -3) {

						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -1;
					}
					else if (layer_ptr[tmp.metal_layer - 1][y][x - 1].sym == -4) {

						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -3;
					}
				}

			}


			else {
				if (y + 1 < max_y + 1)
				{
					if (layer_ptr[tmp.metal_layer - 1][y + 1][x].sym == -3) {

						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -2;
					}

					else if (layer_ptr[tmp.metal_layer - 1][y + 1][x].sym == -4) {

						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -3;
					}
				}
				if (y - 1 >= 0) {
					if (layer_ptr[tmp.metal_layer - 1][y - 1][x].sym == -3) {

						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -2;
					}
					else if (layer_ptr[tmp.metal_layer - 1][y - 1][x].sym == -4) {

						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -3;
					}
				}
			}
		}

		for (int i = serials + 1; i<multi_pin_net.size(); i++) {
			for (int j = 0; j<multi_pin_net[i].connect.size(); j++) {
				pin tmp = multi_pin_net[i].connect[j];
				int y = tmp.y * 2;
				int x = tmp.x * 2;
				if (tmp.metal_layer - 1 == 0 || tmp.metal_layer - 1 == 2) // -------layer
				{

					if (x + 1 < max_x + 1)
					{
						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -3;
					}
					if (x - 1 >= 0) {
						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -3;
					}
				}
				else {

					if (y + 1 < max_y + 1)
					{
						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -3;
					}

					if (y - 1 >= 0) {
						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -3;
					}

				}
			}
		}

		for (int i = 0; i<net_Data.size(); i++) {
			for (int j = 0; j<net_Data[i].connect.size(); j++) {
				pin tmp = net_Data[i].connect[j];
				int y = tmp.y * 2;
				int x = tmp.x * 2;
				if (tmp.metal_layer - 1 == 0 || tmp.metal_layer - 1 == 2) // -------layer
				{

					if (x + 1 < max_x + 1)
					{
						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -3;

					}

					if (x - 1 >= 0) {
						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -3;
					}
				}
				else {

					if (y + 1 < max_y + 1)
					{
						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -3;
					}

					if (y - 1 >= 0) {
						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -3;
					}

				}
			}
		}
	}


	else {

		for (int j = 0; j < net_Data[serials].connect.size(); j++) {
			pin tmp = net_Data[serials].connect[j];
			int y = tmp.y * 2;
			int x = tmp.x * 2;
			if (tmp.metal_layer - 1 == 0 || tmp.metal_layer - 1 == 2) {

				if (x + 1 < max_x + 1)
				{
					if (layer_ptr[tmp.metal_layer - 1][y][x + 1].sym == -3)
						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -1;

					else if (layer_ptr[tmp.metal_layer - 1][y][x + 1].sym == -4)
						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -3;
				}
				if (x - 1 >= 0) {
					if (layer_ptr[tmp.metal_layer - 1][y][x - 1].sym == -3)
						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -1;
					else if (layer_ptr[tmp.metal_layer - 1][y][x - 1].sym == -4)
						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -3;
				}

			}

			else {
				if (y + 1 < max_y + 1)
				{
					if (layer_ptr[tmp.metal_layer - 1][y + 1][x].sym == -3)
						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -2;

					else if (layer_ptr[tmp.metal_layer - 1][y + 1][x].sym == -4)
						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -3;
				}
				if (y - 1 >= 0) {
					if (layer_ptr[tmp.metal_layer - 1][y - 1][x].sym == -3)
						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -2;
					else if (layer_ptr[tmp.metal_layer - 1][y - 1][x].sym == -4)
						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -3;
				}
			}

		}

		for (int i = serials + 1; i<net_Data.size(); i++) {
			for (int j = 0; j<net_Data[i].connect.size(); j++) {
				pin tmp = net_Data[i].connect[j];
				int y = tmp.y * 2;
				int x = tmp.x * 2;

				if (tmp.metal_layer - 1 == 0 || tmp.metal_layer - 1 == 2) // -------layer
				{

					if (x + 1 < max_x + 1)
					{
						layer_ptr[tmp.metal_layer - 1][y][x + 1].sym = -3;

					}

					if (x - 1 >= 0) {
						layer_ptr[tmp.metal_layer - 1][y][x - 1].sym = -3;
					}
				}
				else {

					if (y + 1 < max_y + 1)
					{
						layer_ptr[tmp.metal_layer - 1][y + 1][x].sym = -3;
					}

					if (y - 1 >= 0) {
						layer_ptr[tmp.metal_layer - 1][y - 1][x].sym = -3;
					}
				}
			}
		}
	}
}

void two_pin_routing()
{


	int index = 0;
	while(net_Data[index].critical == 'Y')
		index++;

	vector< vector<int> > order;

	order = rand_order(net_Data.size() - index, 1);

	for(int i=0;i<order[0].size();i++){
		order[0][i]+=index;
	}

	//cout << "critical" << endl;
	for(int i = 0;i < index; i++){
		//


		//cout << net_Data[i].serial_number << endl;
		init_value(net_Data[i].serial_number, 'S', i);

		grid s(net_Data[i].connect[0].metal_layer - 1, (net_Data[i].connect[0].y * 2), (net_Data[i].connect[0].x * 2));

		s.x_point = net_Data[i].connect[0].x;
		s.y_point = net_Data[i].connect[0].y;
		grid t(net_Data[i].connect[1].metal_layer - 1, (net_Data[i].connect[1].y * 2), (net_Data[i].connect[1].x * 2));

		t.x_point = net_Data[i].connect[1].x;
		t.x_point = net_Data[i].connect[1].y;

		layer_ptr[s.i][s.j][s.k].conbine = true;
		layer_ptr[t.i][t.j][t.k].conbine = true;
		int dist = Pioneer(s);
		find_grid(net_Data[i].serial_number,  net_Data[i].serial_number);
	}



	//cout << "Non critical" << endl;
	for (int i = 0; i < net_Data.size() - index; i++) {
		//cout << net_Data[order[0][i]].serial_number << endl;
		init_value(net_Data[order[0][i]].serial_number, 'S',order[0][i]);

		grid s(net_Data[order[0][i]].connect[0].metal_layer - 1, (net_Data[order[0][i]].connect[0].y * 2), (net_Data[order[0][i]].connect[0].x * 2));

		s.x_point = net_Data[order[0][i]].connect[0].x;
		s.y_point = net_Data[order[0][i]].connect[0].y;
		grid t(net_Data[order[0][i]].connect[1].metal_layer - 1, (net_Data[order[0][i]].connect[1].y * 2), (net_Data[order[0][i]].connect[1].x * 2));

		t.x_point = net_Data[order[0][i]].connect[1].x;
		t.x_point = net_Data[order[0][i]].connect[1].y;

		layer_ptr[s.i][s.j][s.k].conbine = true;
		layer_ptr[t.i][t.j][t.k].conbine = true;
		int dist = Pioneer(s);
		find_grid(net_Data[order[0][i]].serial_number,  net_Data[order[0][i]].serial_number);
	}

	//cout << "end of non critical" << endl;


}

vector<vector <int> > rand_order(int n, int times)
{


	int *a = new int[n];
	for (int i = 0; i < n; i++)
		a[i] = i;

	vector< vector<int> > r;
	vector<int> v(a, a + n);
	map<vector<int>, bool> m;

	while (m.size() < times)
	{
		for (int i = n - 1; i > 0; i--)
			swap(v[i], v[rand() % (i + 1)]);
		m[v] = true;
	}

	//cout << "test:" << endl;
	map<vector<int>, bool>::iterator it;
	for (it = m.begin(); it != m.end(); ++it) {
		//cout << "************" << it->first.size() << endl;
		r.push_back(it->first);
	}

	return r;
}


void mul_pin_routing() {

	for (int i = 0; i <multi_pin_net.size(); i++)
	{
		//cout << "net "<<multi_pin_net[i].serial_number<<endl;
		int f = 0;

		vector< vector<int> > order;

		order = rand_order(multi_pin_net[i].edges.size(), 1);
		//cout << order[0].size() << endl;
		//cout << "size: " << multi_pin_net[i].edges.size() << endl;


		/*
		for (int times = 0; times < order[0].size(); times++) {
		cout << order[0][times] << ' ';
		}
		cout << endl;
		*/
		for (int j = 0; j < multi_pin_net[i].edges.size(); j++) {
			//cout << j << "result in: " << order[0][j] << endl;
			int start = multi_pin_net[i].edges[order[0][j]].start_serial;
			int end = multi_pin_net[i].edges[order[0][j]].end_serial;

			grid s(pin_Data[start - 1].metal_layer - 1, (pin_Data[start - 1].y * 2), (pin_Data[start - 1].x * 2));

			s.x_point = pin_Data[start - 1].x;
			s.y_point = pin_Data[start - 1].y;

			grid t(pin_Data[end - 1].metal_layer - 1, (pin_Data[end - 1].y * 2), (pin_Data[end - 1].x * 2));

			t.x_point = pin_Data[end - 1].x;
			t.y_point = pin_Data[end - 1].y;

			init_value(multi_pin_net[i].serial_number, 'M', i);
			if (f == 0) {
				layer_ptr[t.i][t.j][t.k].conbine = true;
				layer_ptr[s.i][s.j][s.k].conbine = true;
				f++;
				int dist1 = Pioneer(s);
			}
			else {
				int dist1 = Pioneer(s);
				init_travelsal();
				int dist2 = Pioneer(t);
			}
			//cout << "end\n" << endl;
		}
		find_grid(multi_pin_net[i].serial_number, multi_pin_net[i].serial_number);

	}
}

int main(int argc, char*argv[])
{
	srand(time(NULL));
	double start, end;
	int contiues = 1;

	if (argc < 5) {
		return false;
	}
	for (int i = 1; i < 4; i++) {
		readFile(argv[i], i - 1);   //"case1_pin.in","case1_net.in","case1_blockage.in"
	}
	max_x *= 2;
	max_y *= 2;
	//bool enhance=true;
	int run = 0;
	double pre_scores = (double)INT_MAX;
	bool average=true;
	start = clock();
	double average_length=0.0;
	double average_via=0.0;
	double average_Detour_critical=0.0;

	while (contiues) {
		global_flag = 1;
		total_wire_length = 0;
		total_via = 0;
		Detour_critical = 0;
		balance_score = 0;
		find_error = false;
		init_circuit();
		mul_pin_routing();
		two_pin_routing();
		vertify();
		//cout << "after vertify" << endl;
		double score = 0.0;
		if (global_flag && !find_error) {
			count_balance_ratio();
			if(average){
				average=false;
				average_length=total_wire_length;
				average_via=total_via;
				average_Detour_critical=Detour_critical;
			}

			//score = (total_wire_length*0.35) + (total_via*0.2) + (Detour_critical*0.35) + ((100 - (balance_score*0.25))*0.1);
			score=(((double)total_wire_length/(double)average_length) * 500  +  ((double)total_via/ (double)average_via )*10 + ((double)Detour_critical/ (double)average_Detour_critical )*500
				+ ((100 - (balance_score*0.25))/100)*10) ;

			if(score<pre_scores){
				cout << "total_wire_length: " << total_wire_length << endl;
				cout << "total_via: " << total_via << endl;
				cout << "Detour_critical: " << Detour_critical << endl;
				cout << "balance_score: " << balance_score * 0.25 << endl;
				cout << fixed << setprecision(2) << "Final score: " << score << endl << endl;
				pre_scores = score;
				copy_net_of_storage.clear();
				for (iter = net_of_storage.begin(); iter != net_of_storage.end(); iter++) {
					copy_net_of_storage.insert(pair<int, vector<clear_net> >(iter->first,iter->second));
				}

				ofstream out(argv[argc - 1], ios::out);
				print_net_of_storage(out, pre_scores,argv[argc - 1]);

			}
		}

		net_of_storage.clear();
		balance_ratio.clear();
		//out.close();
		end = clock();
		double time = double(end - start) / CLOCKS_PER_SEC;
		//cout << "Cost time: " << time << " ms" << endl << endl;
		if (time > 7000) {
			contiues = 0;
		}


	}

    return 0;
}
