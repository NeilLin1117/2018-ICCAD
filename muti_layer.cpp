#include <iostream>
#include <list>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <sstream>
#include <fstream>
using namespace std;

class grid
{

public:
	grid():value(0),visited(false){}

	grid(int i,int j,int k){
		value = 0;
		visited = false;
		this->i = i;
		this->j = j;
		this->k = k;

	}
	//we know that 0,2 is horizontal
	//and 1,3 is vertical
	int i;
	int j;
	int k;
	bool visited;
	int value;
	int fi;
	int fj;
	int fk;
	char sym;

};

class pin {
public:
	pin():x(0),y(0),serial_number(0),metal_layer(0){}
	pin &operator = (const pin &inpin){
		x = inpin.x;
		y = inpin.y;
		serial_number = inpin.serial_number;
		metal_layer = inpin.metal_layer;
	}
	int serial_number;
	int metal_layer;
	int x;
	int y;
};

class net {
public:
	net():priority(0){}
	int serial_number;
	vector <pin> connect;
	char critical;
	int priority;
	bool operator<(const net &rhs) const { return priority < rhs.priority; }


};

class blockage {
public:
	blockage():x1(0),y1(0),x2(0),y2(0){}
	int x1;            //¥ª¤Ux§¤¼Ð
	int y1;            //¥ª¤Uy§¤¼Ð
	int x2;            //¥k¤Wx§¤¼Ð
	int y2;            //¥k¤Wy§¤¼Ð
};

vector <string> string_data; //ÅªÀÉ
vector <pin> pin_Data;
vector <net> net_Data;
vector <blockage> blk_Data;
int max_x = 0;
int max_y = 0;



//gobal decalre
grid ***layer_ptr = new grid **[4];
vector <grid> find_neighbor(grid point);
void traceback(grid s,grid t);
int count_words(string s);
int count_field_num(pin p1,pin p2);
void init_value();


int BFS(grid src,grid tgt)
{
	list<grid>queue;
	layer_ptr[src.i][src.j][src.k].visited = true;
	queue.push_back(src);


	cout << "src: " <<'('<< src.i <<','<< src.j <<','<<src.k <<')'<<endl;

	while(!queue.empty()){
		src = queue.front();
		//cout << src.i << ',' << src.j << ' ' << src.k <<endl;
		queue.pop_front();


		if(src.i == tgt.i && src.j == tgt.j && src.k == tgt.k){
			cout << "find point" << endl;
			return layer_ptr[src.i][src.j][src.k].value;
			break;
		}
		vector <grid> neighbor = find_neighbor(src);


		for(int i=0;i<neighbor.size();i++){
			int li = neighbor[i].i;
			int lj = neighbor[i].j;
			int lk = neighbor[i].k;
			//layer_ptr[li][lj][lk].visited == false
			if(layer_ptr[li][lj][lk].visited == false){

				layer_ptr[li][lj][lk].visited = true;
				layer_ptr[li][lj][lk].value = src.value + 1;
				layer_ptr[li][lj][lk].fi = src.i;
				layer_ptr[li][lj][lk].fj = src.j;
				layer_ptr[li][lj][lk].fk = src.k;
				queue.push_back(layer_ptr[li][lj][lk]);
			}
		}
	}
}

void traceback(grid s,grid t)
{

	int tmpi = t.i;
	int tmpj = t.j;
	int tmpk = t.k;
	char net = layer_ptr[t.i][t.j][t.k].sym;

	while(tmpi != s.i || tmpj != s.j || tmpk != s.k ){

		layer_ptr[t.i][t.j][t.k].sym = net;
		tmpi = layer_ptr[t.i][t.j][t.k].fi;
		tmpj =  layer_ptr[t.i][t.j][t.k].fj;
		tmpk =  layer_ptr[t.i][t.j][t.k].fk;

		t.i = tmpi;
		t.j = tmpj;
		t.k =tmpk;
		cout << "ti : point is(" << t.i << ',' << t.j << ',' << t.k <<')' << endl;

	}
}

vector <grid> find_neighbor(grid point)
{
	vector <grid> neighbor;

	if(point.i % 2 == 0) {// the layer horizontal

		if(point.k+1 < max_x + 1){
			//layer_ptr[point.i][point.j][point.k+1].sym = '-';
			if(layer_ptr[point.i][point.j][point.k+1].sym == '-')
				neighbor.push_back(layer_ptr[point.i][point.j][point.k+1]);

		}
		if(point.k-1 > 0){
			//layer_ptr[point.i][point.j][point.k-1].sym = '-';
			if(layer_ptr[point.i][point.j][point.k-1].sym == '-')
				neighbor.push_back(layer_ptr[point.i][point.j][point.k-1]);
		}

	}
	else{ // the layer is vertical

		if(point.j + 1 < max_y + 1){
			//layer_ptr[point.i][point.j][point.k+1].sym = '|';
			if(layer_ptr[point.i][point.j + 1][point.k].sym == '|')
				neighbor.push_back(layer_ptr[point.i][point.j + 1][point.k]);
		}
		if(point.j - 1 > 0){
			//layer_ptr[point.i][point.j][point.k-1].sym = '|';
			if(layer_ptr[point.i][point.j - 1][point.k].sym == '|')
				neighbor.push_back(layer_ptr[point.i][point.j - 1][point.k]);
		}

	}
	if(point.i + 1 < 4){
		if(layer_ptr[point.i + 1][point.j][point.k].sym != 'x')
			neighbor.push_back(layer_ptr[point.i + 1][point.j][point.k]);
	}
	if(point.i - 1 >= 0){
		if(layer_ptr[point.i - 1][point.j][point.k].sym != 'x')
			neighbor.push_back(layer_ptr[point.i - 1][point.j][point.k]);
	}

	if(neighbor.empty()){
		cout << "lock no route!!" << endl;
		exit(1);
	}

	return neighbor;
}

void Parser(int choice)
{
	if (choice == 0 || choice == 2) {

		for (int i = 0; i<string_data.size(); i++) {

			stringstream ss(string_data[i]);

			if (choice == 0) {


				pin a;
				ss >> a.serial_number;
				ss >> a.metal_layer;
				ss >> a.x;
				ss >> a.y;

				if(a.x > max_x)
					max_x = a.x;
				if(a.y > max_y)
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

			stringstream ss(string_data[i]);
			net a;

			ss >> a.serial_number;
			int tmp;

			for(int i=0;i < str_size - 1;i++) {
				ss >> tmp;
				a.connect.push_back(pin_Data[tmp]);
			}

			ss >> a.critical;


			for (int j = 0; j < a.connect.size() - 1 ; j++) {

				int sub_x = abs(a.connect[j].x - a.connect[j+1].x);
				int sub_y = abs(a.connect[j].y - a.connect[j+1].y);
				a.priority += sub_x + sub_y;
			}


			net_Data.push_back(a);

		}

	}

}


int count_words(string s)
{
	int words = 0;
	for(int i = 0;i < s.length();i++){
		if(s[i] == ' '){
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
	}   // end if
	//char s[200];
	string_data.clear();
	string tmp;
	while (getline(infile,tmp)) {

		//getline(infile,tmp); //read name of document
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


void print()
{


	for(int i=0;i<4;i++){
		for(int j=max_y ;j>=0;j--){
			for(int k=0;k<max_x + 1;k++){
				cout   << layer_ptr[i][j][k].sym ;
			}
			cout << endl;
		}
		cout << 'l' << endl;
	}

}



void init_circuit()
{

	//for(int i=0;i<max_y)


	for(int i=0;i < 4;i++){
		layer_ptr[i] = new grid* [max_y + 1];
		for(int j=0;j < max_y + 1;j++)
			layer_ptr[i][j] = new grid[max_x + 1]();
	}

	for(int i=0;i<4;i++){

		if(i % 2 == 0){
			for(int j=0;j<max_y + 1;j++){
				for(int k=0;k<max_x + 1;k++){
					layer_ptr[i][j][k].sym = '-';
					layer_ptr[i][j][k].i = i;
					layer_ptr[i][j][k].j = j;
					layer_ptr[i][j][k].k = k;
				}
			}
		}
		else{
			for(int j=0;j<max_y + 1;j++){
				for(int k=0;k<max_x + 1;k++){
					layer_ptr[i][j][k].sym = '|';
					layer_ptr[i][j][k].i = i;
					layer_ptr[i][j][k].j = j;
					layer_ptr[i][j][k].k = k;
				}
			}
		}

	}


	for(int i=0;i<blk_Data.size();i++)
		for(int j=blk_Data[i].y1;j < blk_Data[i].y2;j++)
			for(int k = blk_Data[i].x1;k < blk_Data[i].x2;k++)
				for(int i=0;i<4;i++)
					layer_ptr[i][j][k].sym = 'x';


	for(int i=0;i<net_Data.size();i++)
		for(int j=0;j<net_Data[i].connect.size();j++){
			pin tmp = net_Data[i].connect[j];


			layer_ptr[tmp.metal_layer - 1][tmp.y][tmp.x].sym = net_Data[i].serial_number + '0';

		}

}

void init_value()
{
	for(int i=0;i<4;i++){
		for(int j=max_y ;j>=0;j--){
			for(int k=0;k<max_x + 1;k++){
				layer_ptr[i][j][k].value = 0 ;
				layer_ptr[i][j][k].visited = false;
			}
			//cout << endl;
		}
		//cout << 'l' << endl;
	}
}

void routing(int i)
{
	for(int j = 0; j < net_Data[i].connect.size() - 1;j++){
		grid s(net_Data[i].connect[j].metal_layer - 1,net_Data[i].connect[j].y,net_Data[i].connect[j].x);
		grid t(net_Data[i].connect[j+1].metal_layer - 1,net_Data[i].connect[j+1].y,net_Data[i].connect[j+1].x);

		int dist = BFS(s,t);
		traceback(s,t);
		//print();
		init_value();
	}
}

int main(int argc,char*argv[])
//int main()
{

	if (argc < 4) {
		return false;
	}


	pin tmp;
	pin_Data.push_back(tmp);

	for (int i = 1; i<4; i++) {
		readFile(argv[i], i-1);   //"case1_pin.in","case1_net.in","case1_blockage.in"
	}
	//max_y = 150;
	//max_x = 210;
	init_circuit();
	cout << "the max x is " << max_x <<endl;
	cout << "the max y is " << max_y << endl;

	for(int i=0;i<net_Data.size();i++)
		routing(i);

	print();



}
