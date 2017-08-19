#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
using namespace std;

const double TL = 1.0;
const double TR = 2.0;

enum Rate_Func {
	rf1,
	rf2
};

Rate_Func rate_func = rf1;

typedef list<double> List;

struct interaction
{
	double time;
	int location;
	int level;
	interaction* left;
	interaction* right;
};

struct Profile {
	//double value = 0;
	List values;
	double pre_time = 0;
	double summation = 0;
};

inline double rate_function(double x, double y) {
	if(x >= 0 && y >= 0)
		return (rate_func == rf1) ? sqrt(x*y/(x+y)) : sqrt(x+y);
	else {
		cout<<"error, sqrt of negative number! "<<endl;
		return 0;
	}
}

void print_v(double* Array, int size)
{
	for(int i = 0; i < size; i++)
		cout<< Array[i] << "  ";
	cout<<endl;
}

int find_min(interaction* Link)
{
		double tmp = Link->time;
		interaction* pt = Link;
		interaction* tmp_pt = Link;
		while(tmp_pt != NULL)
		{
			if(tmp_pt->time < tmp)
			{
				tmp = tmp_pt->time;
				pt = tmp_pt;
			}
			tmp_pt = tmp_pt->right;
		}

		return pt->location;

}

void remove(interaction** Link, interaction* pt)
//remove pt from link but keep pt for future use
{
	if( *Link == NULL || pt == NULL)
		return;
	else if( pt->left == NULL && pt->right == NULL)
	{
		*Link = NULL;
		return;
	}
	else
	{
		if((*Link) == pt )
			*Link = pt->right;
		if( pt->right != NULL)
			pt->right->left = pt->left;
		if( pt->left != NULL)
			pt->left->right = pt->right;
	}
	pt->left = NULL;
	pt->right = NULL;
}

void push_front(interaction** Link, interaction* pt)
//push the interaction pointed by pt into the front of the lise
{
	pt->left = NULL;
	if(*Link == NULL)
	{
		pt->right = NULL;
		*Link = pt;
	}
	else
	{
		pt->right = *Link;
		(*Link)->left = pt;
		*Link = pt;
	}
}

void print_list(interaction* Link)
{
	interaction* tmp = Link;
	while(tmp!= NULL)
	{
		cout<<" location: "<< tmp->location <<" time: " << tmp->time << "  " ;
		tmp = tmp->right;
	}
	cout<<endl;
}


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, const int N, const double small_tau, const int ratio, const int Step)
//distribute clock times of a big step into vectors that represent small steps. If the clock time is bigger than a big tau, then it is arranged in the right location
{
	int tmp;
	for(int i = 0; i < N; i++)
	{
		if(time_array[i].time > (Step + 1)*ratio*small_tau)
		{
			tmp = ratio;
		}
		else
		{
			tmp = floor((time_array[i].time - ratio*small_tau*Step)/small_tau);
			if(tmp < 0) {
				cout<<"Error! tmp = "<<tmp<<" when i = "<<i<<" small tau = "<<small_tau<<endl;
				cout<<"time difference is "<<time_array[i].time - ratio*small_tau*Step<<endl;
				exit(1);
			}
		}
		
		time_array[i].level = tmp;
		push_front(&clock_time_in_step[tmp], &time_array[i] );

	}
}



void move_interaction(interaction** &clock_time_in_step, interaction* pt, const double small_tau, const int ratio, const int Step, const double new_time)
//move the interaction pointed by *pt from old bucket to new bucket
//n_move1: move without relinking pointers
//n_move2: move with relinking pointers
{
	//double old_time = pt->time;
	int old_level, new_level;
	pt->time = new_time;
	old_level = pt->level;

	if ( new_time > (Step + 1)*ratio*small_tau )
	{
		new_level = ratio;
	}
	else
	{
		new_level = floor( (new_time - Step*ratio*small_tau)/small_tau );
		if( new_level < 0 || new_level >= ratio)
		{
			cout<<"Error!!"<<endl;
			cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
			cout<<"error time is "<<new_time<<endl;
			exit(1);
		}
	}
	//    cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
	if(old_level != new_level ) {
		//pt->time = new_time;
		pt->level = new_level;
		remove(&(clock_time_in_step[old_level]), pt);
		push_front(&(clock_time_in_step[new_level]), pt);
	}
}



void update(interaction** &clock_time_in_step, const int level, const int N, const double small_tau, const int ratio, const int Step, interaction* time_array, double *energy_array, trng::uniform01_dist<> &u, trng::yarn2 &r, int &count, Profile &X, Profile &Y, Profile &X_times_Y, const int mid_point)
{
	double next_time = (Step*ratio + level + 1)*small_tau;
	int min_loc = find_min(clock_time_in_step[level]);
	interaction* pt = &time_array[min_loc];
	double current_time = pt->time;
	while(clock_time_in_step[level] != NULL)
	{
		count++;
		//Step 1: update min interaction and energy
		double total_energy = energy_array[min_loc] + energy_array[min_loc + 1];
		double tmp_double;
		double rnd;
		do{
			rnd = u(r);
		} while(rnd < 1e-15 || rnd > 1 -  1e-15);
		tmp_double = - log(rnd)/rate_function(energy_array[min_loc], energy_array[min_loc + 1]);
		if(tmp_double > 1e10 || tmp_double < 1e-15) {
			cout<<"tmp_double ="<<tmp_double<<endl;
			exit(1);
		}
		double old_e_left = energy_array[min_loc];
		double old_e_right = energy_array[min_loc + 1];
		double tmp_rnd_uni;
		do{
			tmp_rnd_uni = u(r);
		} while(tmp_rnd_uni < 1e-15 || tmp_rnd_uni > 1 -  1e-15);

		if(min_loc == 0)
		{
			do{
				rnd = u(r);
			} while(rnd < 1e-15 || rnd > 1 -  1e-15);
			total_energy = old_e_right  -log(rnd)*TL;
		}
		if(min_loc == N)
		{
			do{
				rnd = u(r);
			} while(rnd < 1e-15 || rnd > 1 -  1e-15);
			total_energy = old_e_left  -log(rnd)*TR;
		}
		
		if(min_loc != 0)
		{
			energy_array[min_loc] = tmp_rnd_uni*total_energy;
//            E_avg[min_loc - 1] += (current_time - last_update[min_loc - 1])*old_e_left;
//            last_update[min_loc - 1] = current_time;
		}
		if(min_loc != N)
		{
			energy_array[min_loc + 1] = (1 - tmp_rnd_uni)*total_energy;//update energy
//            E_avg[min_loc] += (current_time - last_update[min_loc])*old_e_right;
//            last_update[min_loc] = current_time;
		}
		
		//do the integration
		//Both energy[int(N/2)] and energy[int(N/2) + 1] are changing
		if(min_loc == mid_point) {
			
			X.summation = (current_time - X.pre_time) * X.values.back();
			X.pre_time = current_time;
			X.values.push_back(energy_array[min_loc]);
			
			Y.summation = (current_time - Y.pre_time) * Y.values.back();
			Y.pre_time = current_time;
			Y.values.push_back(energy_array[min_loc + 1]);
			
			X_times_Y.summation = (current_time - X_times_Y.pre_time)*(X_times_Y.values.back());
			X_times_Y.pre_time = current_time;
			X_times_Y.values.push_back(energy_array[min_loc] * energy_array[min_loc + 1]);
		}
		
		//Only energy[int(N/2)] is changing
		else if (min_loc + 1 == mid_point) {
			
			X.summation = (current_time - X.pre_time) * X.values.back();
			X.pre_time = current_time;
			X.values.push_back(energy_array[min_loc]);
			
			X_times_Y.summation = (current_time - X_times_Y.pre_time)*(X_times_Y.values.back());
			X_times_Y.pre_time = current_time;
			X_times_Y.values.push_back(energy_array[min_loc] * energy_array[min_loc + 1]);
		}
		
		//Only energy[int(N/2) + 1] is changing
		else if (min_loc == mid_point + 1) {
			
			Y.summation = (current_time - Y.pre_time) * Y.values.back();
			Y.pre_time = current_time;
			Y.values.push_back(energy_array[min_loc + 1]);
			
			X_times_Y.summation = (current_time - X_times_Y.pre_time)*(X_times_Y.values.back());
			X_times_Y.pre_time = current_time;
			X_times_Y.values.push_back(energy_array[min_loc] * energy_array[min_loc + 1]);
			
		}
		
		move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);

		//
		//Step 2: update other interactions
		if(min_loc != 0)
		{
			pt = &time_array[min_loc - 1];
			tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc - 1], old_e_left)/rate_function(energy_array[min_loc - 1], energy_array[min_loc]) + current_time;
			if(tmp_double < current_time)
			{
				cout<<"Error with time at "<<min_loc<<" with diff = "<<tmp_double-current_time<<endl;
				cout<<"energy :"<<energy_array[min_loc - 1]<<" "<<energy_array[min_loc]<<endl;
				exit(1);
				
			}
			move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
		}
		if(min_loc != N)
		{
			pt = &time_array[min_loc + 1];
			tmp_double = (pt->time - current_time)*rate_function(energy_array[min_loc + 2], old_e_right)/rate_function(energy_array[min_loc + 2], energy_array[min_loc + 1]) + current_time;
			if(tmp_double < current_time)
			{
				cout<<"Error with time at "<<min_loc<<" with diff = "<<tmp_double-current_time<<endl;
				cout<<"energy :"<<energy_array[min_loc + 1]<<" "<<energy_array[min_loc + 2]<<endl;
				exit(1);
			}
			move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);

		}

		//Step 3: update current time
		if(clock_time_in_step[level] != NULL)
		{
			min_loc = find_min(clock_time_in_step[level]);
			pt = &time_array[min_loc];
			current_time = pt->time;
			if(current_time > next_time + 1e-9 || current_time < next_time - small_tau)
            {
        			cout.precision(17);
                cout<<"Error ! current_time = "<<current_time<<" next_time  = "<<next_time<<" diff = "<<current_time - next_time<<endl;
                cout<<"min_loc = "<<min_loc<<" clock_time_in_step[level] = "<<clock_time_in_step[level]->time<<endl;
                cout<<" level = "<<level<<" count = "<<count<<endl;
				for(int i = 0; i < ratio + 1; i++) {
				    cout<<"level "<<i<<" ";
				    print_list(clock_time_in_step[i]);
				}
				print_v(energy_array, N+2);
				for(int i = 0; i < N+1; i++)
				  cout<<time_array[i].time<<" ";
				cout<<endl;
				exit(1);
				//break;
            }
		}


	}

}

double variance(List data) {
	double n = 0;
	double Sum = 0;
	double SumSq = 0;
	for(double item : data) {
		n += 1;
		Sum += item;
		SumSq += (item * item);
	}
	return (SumSq - (Sum * Sum)/n)/(n-1);
}

int main(int argc, char** argv)
{
	struct timeval t1, t2;
	ofstream myfile;
	myfile.open("HL_KMP.txt", ios_base::app);

	int N = 100;
	if(argc > 1)
	{
		N = strtol(argv[1], NULL,10 );
	}
	double big_tau = 0.2;//big time step of tau leaping
	const int ratio = int(N/10);//ratio of big step and small step
	double small_tau = big_tau/double(ratio);//small time step
	double rnd;

	double* energy_array = new double[N+2];
	double *E_avg = new double[N];
	double* last_update = new double[N];
	for(int i = 0; i <N; i++)
	{
		E_avg[i] = 0;
		last_update[i] = 0;
	}
	trng::yarn2 r;
    r.seed(time(NULL));
	trng::uniform01_dist<> u;
	energy_array[0] = TL;
	energy_array[N+1] = TR;
	for(int n = 1; n < N+1; n++)
		energy_array[n] = 1;
	interaction* time_array = new interaction[N+1];
	for(int n = 0; n < N+1; n++)
	{
		do{
			rnd = u(r);
		} while(rnd < 1e-15 || rnd > 1 -  1e-15);
		time_array[n].time = -log(rnd)/rate_function(energy_array[n], energy_array[n+1]);
		time_array[n].location = n;
		time_array[n].level = -1;
		time_array[n].left = NULL;
		time_array[n].right = NULL;
	}
	
	//X := energy_array[int(N/2)], Y := energy_array[int(N/2) + 1];
	//cov(X,Y) = E[XY]-E[X]E[Y]
	Profile X;
	Profile Y;
	Profile X_times_Y;
	
	X.values.push_back(1);
	Y.values.push_back(1);
	X_times_Y.values.push_back(X.values.back() * Y.values.back());
	
	int count = 0;

	interaction** clock_time_in_step = new interaction*[ratio + 1];//each element in the array is the head of a list
	for(int i = 0; i < ratio + 1; i++)
	{
		clock_time_in_step[i] = NULL;
	}

	gettimeofday(&t1,NULL);

	int Step = 50;

	for(int out_n = 0; out_n < Step; out_n++)
	{
		big_step_distribute(clock_time_in_step,time_array,N+1,small_tau,ratio,out_n);

		for(int in_n = 0; in_n < ratio; in_n++)
		{

			if(clock_time_in_step[in_n]!= NULL)
			{
				update(clock_time_in_step, in_n, N, small_tau, ratio, out_n, time_array, energy_array, u, r, count, X, Y, X_times_Y, int(N/2));
			}
		}

		clock_time_in_step[ratio] = NULL;

	}
	
	//(Step*big_tau) is total time. it will be used to find the probability for expectation
	//∑(energy)*(time_lapsed)/total_time, where (time_lapsed)/total_time is the probability
	
	//Expectation Value Of XY, X, Y (E[XY], E[X], E[Y]);
	cout<<"E[XY] = "<<X_times_Y.summation/(Step*big_tau)<<endl;
	cout<<"E[X] = "<<X.summation/(Step*big_tau)<<endl;
	cout<<"E[Y] = "<<Y.summation/(Step*big_tau)<<endl;
	double cov = X_times_Y.summation/(Step*big_tau) - (X.summation/(Step*big_tau) * Y.summation/(Step*big_tau));
	cout<<"cov(X,Y) = E[XY] - E[X]E[Y] = "<<cov<<endl;
	double varX = variance(X.values);
	double varY = variance(Y.values);
	cout<<"Var(X) = "<<varX<<endl;
	cout<<"Var(Y) = "<<varY<<endl;
	double correlation = cov/(sqrt(varX)*sqrt(varY));
	cout<<"Correlation = "<<correlation<<endl;
//	for(int i = 0 ; i < N + 2 ; i++) {
//		cout<<summation[i]/(Step*big_tau)<<" ";
//	}
//	cout<<endl;

	gettimeofday(&t2, NULL);
	delete[] energy_array;
	delete[] E_avg;
	delete[] time_array;
	delete[] clock_time_in_step;

	double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
					t2.tv_usec - t1.tv_usec) / 1.e6;

//    cout << "total CPU time = " << delta <<endl;
	cout<<" N = "<<N <<endl;
	//cout<<"seconds per million event is "<< 1000000*delta/double(count)<<endl;
	//myfile<<" N = "<<N <<endl;
	myfile<< 1000000*delta/double(count)<<"  ";
	myfile.close();

}
