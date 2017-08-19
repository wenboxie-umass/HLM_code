#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <list>
#include <vector>
#include <map>
#include <omp.h>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <array>
using namespace std;

typedef pair<int, int> Pair;
typedef vector<Pair> P_list;
typedef vector<int> List;
typedef vector<double> Double_list;
typedef List::iterator iter;
typedef map<Pair, bool> Dictionary;

double TL = 1.0;
double TR = 2.0;

/*
horizontal:
        |-------------|
        |             |
        |             |
  |------------c------------|
  |     |             |     |
  |     c             c     |
  |     |             |     |
  |-----------c_i-----------|
  |     |             |     |
  |     c             c     |
  |     |             |     |
  |------------c------------|
        |             |
        |             |
        |-------------|

vertical:
        |------|------|
        |      |      |
  |--------c-------c--------|
  |     |      |      |     |
  |     c     c_i     c     |
  |     |      |      |     |
  |--------c-------c--------|
        |      |      |
        |------|------|
*/

enum Rate_Func {
  rf1,
  rf2
};

enum Direction {
    horizontal,
    vertical
 };

struct interaction
{
    double time;
    int location; // Store the index in 1D time_array
    int level;
    Pair coordinate;
    List adjacent_clocks;
    Direction direction; //indicates position of clock. it can be either on vertical line or on horizontal line
    interaction* left;
    interaction* right;
};

struct site {
    double previous_time;
    double previous_energy;
    Pair coordinate;
    Double_list record;
    double summation;
};

Rate_Func rate_func = rf1;

//P_list find_coordinate_of_site_in_energy_array(interaction info, int rows, int cols);

int find_index(int i, int j, int N, int M, Direction direction)
{
	if(direction == vertical)
	{
		return (int)((M + 1) * (int)(i / 2) + M * (int)(i / 2) + (j));
	}

	else
	{
		return (int)((M + 1) * (int)(i / 2 + 1) + M * (int)(i / 2) + (j));
	}
}

inline double rate_function(double x, double y) {

    if(x >= 0 && y >= 0)
        return (rate_func == rf1) ? sqrt(x*y/(x+y)) : sqrt(x+y);
    else {
        cout<<"error, sqrt of negative number! "<<endl;
        return 0;
    }
}

//find_adjacent_clocks will return all possible clocks which are adjacent to clock at (i,j).
List find_adjacent_clocks(int i, int j, Direction direction, int N, int M, int current_index) {
    List lst;
    switch (direction) {
      case vertical:
          if(i != 0) {
              if(j != 0) {
                  //we have a clock at (i - 1, j - 1);
                  lst.push_back(find_index(i - 1, j - 1, N, M, horizontal));
              }

              if(j != M) {
                  //we have a clock at (i - 1, j)
                  lst.push_back(find_index(i - 1, j, N, M, horizontal));
              }
          }

          if(i != 2 * N - 2) {
              if(j != 0) {
                  //we have a clock at (i + 1, j - 1);
                  lst.push_back(find_index(i + 1, j - 1, N, M, horizontal));
              }

              if(j != M) {
                  //we have a clock at (i + 1, j)
                  //cout<<i + 1<<","<<j<<", index: "<<find_index(i + 1, 0, N, M, horizontal)<<endl;
                  lst.push_back(find_index(i + 1, j, N, M, horizontal));
              }
          }

          if(j != 0) {
              lst.push_back(find_index(i, j - 1, N, M, vertical));
          }

          if(j != M) {
              lst.push_back(find_index(i, j + 1, N, M, vertical));
          }
          break;
      case horizontal:
          if(i != 1) {
              lst.push_back(find_index(i - 2, j, N, M, horizontal));
          }

          if(i != 2 * N - 3) {
              lst.push_back(find_index(i + 2, j, N, M, horizontal));
          }

          lst.push_back(find_index(i - 1, j, N, M, vertical));
          lst.push_back(find_index(i - 1, j + 1, N, M, vertical));
          lst.push_back(find_index(i + 1, j, N, M, vertical));
          lst.push_back(find_index(i + 1, j + 1, N, M, vertical));
          break;
    }
    //lst.push_back(10);
    return lst;
}

double energy_summation(double **energy_array, int x, int y, Direction direction) {
    double ret = 10.0;
    switch (direction) {
      case vertical:
        ret = rate_function(energy_array[x / 2][y],energy_array[x / 2][y + 1]);
        break;
      case horizontal:
        ret = rate_function(energy_array[int(x / 2)][y + 1], energy_array[int(x / 2) + 1][y + 1]);
        break;
    }
    return ret;
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
//push the interaction pointed by pt into the front of the list
{
    pt->left = NULL;
    if(*Link == NULL)
    {
        pt->right = NULL;
        *Link = pt;
    }
    else
    {
        //cout<<"ok";
        pt->right = *Link;
        //if((*Link)->left == NULL) cout<<"NULL"<<endl;
        (*Link)->left = pt;
         *Link = pt;
    }
}

void print_list(interaction* Link)
{
    interaction* tmp = Link;
    while(tmp!= NULL)
    {
        //cout<<" location: "<< tmp->location <<" time: " << tmp->time << "  " ;
        tmp = tmp->right;
    }
    cout<<endl;
}


void big_step_distribute(interaction** &clock_time_in_step, interaction* time_array, const int N, const double small_tau, const int ratio, const int Step)
//distribute clock times of a big step into vectors that represent small steps. If the clock time is bigger than a big tau, then it is arranged in the right location
{
    for(int i = 0; i < N; i++)
    {
        int tmp;
        if(time_array[i].time > (Step + 1)*ratio*small_tau)
        {
            tmp = ratio;
        }
        else
        {
            tmp = floor((time_array[i].time - ratio*small_tau*Step)/small_tau);
            if(tmp < 0)
            {
              cout<<"Error! tmp = "<<tmp<<" when i = "<<i<<" small tau = "<<small_tau<<endl;
              cout<<"time difference is "<<time_array[i].time - ratio*small_tau*Step<<endl;
              exit(0);
            }
        }

      //  if(tmp < 0 || tmp > ratio) {
      //      tmp = ratio;
      //  }

      //cout<<"Time: "<<time_array[i].time<<" Tmp: "<<tmp<<endl;
        time_array[i].level = tmp;
        push_front(&clock_time_in_step[tmp], &time_array[i] );
        //cout<<tmp<<endl;


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
    old_level = pt -> level;

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
        }
    }
    //    cout<<"start to move "<< pt->location << " from " << old_level << " to " << new_level<<endl;
    if(old_level == new_level )
    {
        pt->time = new_time;
    }
    else
    {
        pt->level = new_level;
        remove(&(clock_time_in_step[old_level]), pt);
        push_front(&(clock_time_in_step[new_level]), pt);
    }
}



void update(interaction** &clock_time_in_step, const int level, const int N, const int M, const double small_tau, const int ratio, const int Step, interaction* time_array, double **energy_array, trng::uniform01_dist<> &u, trng::yarn2 &r, int *count, double &energy_integration)
//update clock_time_in_step[level]
{
    //cout<<"Start!"<<endl;
    //double next_time = (Step*ratio + level + 1)*small_tau;
    int min_loc = find_min(clock_time_in_step[level]); //find_min should return index of clock in time_array
    //cout<<"OK, Min_Loc: "<<min_loc<<endl;
    interaction* pt = &time_array[min_loc];
    double current_time = pt->time;
//    cout<<"at level " << level << endl;
    Direction direction_min_loc = pt -> direction;
    Pair coordinate_min_loc = pt -> coordinate; // Find the coordinate of min_loc in 2D energy_array;
    int x = coordinate_min_loc.first / 2; //Compute the x coordinate of cell that the min_loc is corresponding to.
    int y = coordinate_min_loc.second;//(direction_min_loc == vertical) ? coordinate_min_loc.second : coordinate_min_loc.second + 1; //Compute the y coordinate of cell......
    List* adjacent_clocks_min_loc = &(pt -> adjacent_clocks);
    Dictionary old_e_value_classifier; // To determine old_e_value I should use to compute the tmp_double
    Pair left = (direction_min_loc == vertical) ? make_pair(x, y) : make_pair(x, y + 1);
    Pair right = (direction_min_loc == vertical) ? make_pair(x, y + 1) : make_pair(x + 1, y + 1);
    old_e_value_classifier[left] = true;
    old_e_value_classifier[right] = true;

    Pair tmp_loc;
    Pair tmp_left, tmp_right;
    int tmp_x = -1;
    int tmp_y = -1;
    Direction dir = vertical;
    double total_energy, tmp_double, old_e_left, old_e_right, tmp_rnd_uni, tmp_rate;

    double pre_energy = 0.0;
    double curr_energy = 0.0;

    double rnd;

    while(true)
    {
        (*count)++;
        total_energy =  (direction_min_loc == vertical) ? energy_array[x][y]+ energy_array[x][y + 1] : energy_array[x][y + 1] + energy_array[x + 1][y + 1];
        tmp_rate = (direction_min_loc == vertical) ? rate_function(energy_array[x][y],energy_array[x][y + 1]) : rate_function(energy_array[x][y + 1] , energy_array[x + 1][y + 1]); //If the direction is horizontal, then we are going to sum the energy of cell at index (x, y) and (x, y + 1). Otherwise, we have to compute the total energy of cell at index (x , y + 1) and (x + 1, y + 1);
        do{
          rnd = u(r);
        }while(rnd < 1e-15 || rnd > 1 -  1e-15);
        tmp_double = -log(rnd)/tmp_rate;

        old_e_left = (direction_min_loc == vertical) ? energy_array[x][y] : energy_array[x][y + 1];
        old_e_right = (direction_min_loc == vertical) ? energy_array[x][y + 1] : energy_array[x + 1][y + 1];
        Pair old_e_left_coordinate, old_e_right_coordinate;
        old_e_left_coordinate = (direction_min_loc == vertical) ? make_pair(x, y) : make_pair(x, y + 1);
        old_e_right_coordinate = (direction_min_loc == vertical) ? make_pair(x, y + 1) : make_pair(x + 1, y + 1);
        do{
            tmp_rnd_uni = u(r);
        }while(tmp_rnd_uni < 1e-15 || tmp_rnd_uni > 1 -  1e-15);

        //store previous energy for energy flux
        if(direction_min_loc == vertical) {
            if(y == 0) {
                pre_energy = energy_array[old_e_right_coordinate.first][old_e_right_coordinate.second];
            }
            else {
                pre_energy = energy_array[old_e_left_coordinate.first][old_e_left_coordinate.second];
            }
        }

        if(y == 0 && direction_min_loc == vertical)
        {
            do{
                rnd = u(r);
            }while(rnd < 1e-15 || rnd > 1 -  1e-15);
            total_energy = old_e_right  -log(rnd)*TL;
        }
        if(y == M && direction_min_loc == vertical)
        {
            do{
              rnd = u(r);
            }while(rnd < 1e-15 || rnd > 1 -  1e-15);
            total_energy = old_e_left  -log(rnd)*TR;
        }

        if(direction_min_loc == vertical) {
            if(y != 0) {
                energy_array[x][y] = tmp_rnd_uni*total_energy;
            }

            if(y != M) {
                energy_array[x][y + 1] = (1 - tmp_rnd_uni)*total_energy;//update energy
            }
        }
        else {
            energy_array[x][y + 1] = tmp_rnd_uni*total_energy;
            energy_array[x + 1][y + 1] = (1 - tmp_rnd_uni)*total_energy;//update energy
        }

        //get current enegy and do the integration
        if(direction_min_loc == vertical) {
          if(y == 0) {
              curr_energy = energy_array[old_e_right_coordinate.first][old_e_right_coordinate.second];
          }
          else {
              curr_energy = energy_array[old_e_left_coordinate.first][old_e_left_coordinate.second];
          }

          energy_integration += (1 - 2*int(y == 0))*(curr_energy - pre_energy);
        }

        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step,  current_time + tmp_double);

        //Step 2: update other interactions
        for(iter i = adjacent_clocks_min_loc->begin() ; i != adjacent_clocks_min_loc->end() ; i++) {
            pt = &time_array[(*i)];
            dir = pt -> direction;
            tmp_loc = pt -> coordinate;
            tmp_x = tmp_loc.first / 2;
            tmp_y = tmp_loc.second;
            if(dir == vertical) {
                //cout<<"Vertical"<<endl;
                tmp_left = make_pair(tmp_x, tmp_y);
                tmp_right = make_pair(tmp_x, tmp_y + 1);
                if(old_e_value_classifier[tmp_left] && direction_min_loc == vertical) {
		  tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_right.first][tmp_right.second],old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                    move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                }

                if(old_e_value_classifier[tmp_right] && direction_min_loc == vertical) {
		  tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_left.first][tmp_left.second], old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                    move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                }

                if(old_e_value_classifier[tmp_left] && direction_min_loc == horizontal) {
                    if(tmp_left == left) {
		      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_right.first][tmp_right.second], old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                    else if(tmp_left == right){
		      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_right.first][tmp_right.second], old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                }

                if(old_e_value_classifier[tmp_right] && direction_min_loc == horizontal) {
                    if(tmp_right == left) {
                      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_left.first][tmp_left.second], old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                      move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                    else if(tmp_right == right){
                      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_left.first][tmp_left.second] , old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                      move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                }
            }
            else {
                tmp_left = make_pair(tmp_x, tmp_y + 1);
                tmp_right = make_pair(tmp_x + 1, tmp_y + 1);
                if(old_e_value_classifier[tmp_right] && direction_min_loc == vertical) {
                    if(tmp_right == left) {
		      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_left.first][tmp_left.second] , old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                    else if(tmp_right == right) {
		      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_left.first][tmp_left.second] , old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                        move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                }

                if(old_e_value_classifier[tmp_left] && direction_min_loc == vertical) {
                    if(tmp_left == left) {
                      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_right.first][tmp_right.second] , old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                      move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                    else if(tmp_left == right) {
                      tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_right.first][tmp_right.second] , old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                      move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                    }
                }

                if(old_e_value_classifier[tmp_right] && direction_min_loc == horizontal) {
		  tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_left.first][tmp_left.second] , old_e_left)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                    move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                }

                if(old_e_value_classifier[tmp_left] && direction_min_loc == horizontal) {
		  tmp_double = (pt->time - current_time)*rate_function(energy_array[tmp_right.first][tmp_right.second] , old_e_right)/energy_summation(energy_array, tmp_loc.first, tmp_loc.second, dir) + current_time;
                    move_interaction(clock_time_in_step, pt,small_tau,ratio, Step, tmp_double);
                }
            }
        }

        //Step 3: update current time
        if(clock_time_in_step[level] != NULL)
        {
            min_loc = find_min(clock_time_in_step[level]);
            pt = &time_array[min_loc];

            current_time = pt->time;

            direction_min_loc = pt -> direction;
            coordinate_min_loc = pt -> coordinate; // Find the coordinate of min_loc in 2D energy_array;
            x = coordinate_min_loc.first / 2; //Compute the x coordinate of cell that the min_loc is corresponding to.
            y = coordinate_min_loc.second;//(direction_min_loc == vertical) ? coordinate_min_loc.second : coordinate_min_loc.second + 1; //Compute the y coordinate of cell......
            adjacent_clocks_min_loc = &(pt -> adjacent_clocks);

            old_e_value_classifier.clear();

            left = (direction_min_loc == vertical) ? make_pair(x, y) : make_pair(x, y + 1);
            right = (direction_min_loc == vertical) ? make_pair(x, y + 1) : make_pair(x + 1, y + 1);
            old_e_value_classifier[left] = true;
            old_e_value_classifier[right] = true;
            //printf("Current Count: %d\n", *count);
        }
        else
        {
            //current_time = next_time + 1;
            //printf("Break!\n");
            break;
        }


    }

}

int main(int argc, char** argv)
{
    struct timeval t1, t2;
    ofstream myfile;
    ofstream data;
    string file_name;

    //Set M = 5 always
    int N = 4, M = 4;
    int N_thread = 8;
    int Step = 100000;
    int type_of_rate_func = 1;
    int current_line_number = -1;
    int current_attempt = -1;
    int number_of_attempts = -1;
    if(argc == 12)
    {
        N = strtol(argv[1], NULL, 0);
        M = strtol(argv[2], NULL, 0);
        N_thread = strtol(argv[3], NULL, 0);
        Step = strtol(argv[4], NULL, 0);
        type_of_rate_func = strtol(argv[5], NULL, 0);
        TL = stod(argv[6]);
        TR = stod(argv[7]);
        file_name = string(argv[8]);
        current_line_number = strtol(argv[9], NULL, 0);
        current_attempt = strtol(argv[10], NULL, 0);
        number_of_attempts = strtol(argv[11], NULL, 0);
    }

    if(type_of_rate_func == 1) {
        rate_func = rf1;
    }
    else {
        rate_func = rf2;
    }

    data.open(file_name, ios_base::app);

    const double big_tau = 0.2;//big time step of tau leaping
    double* energy_integration = new double[N_thread];

    //Test Begin Here

    const int size_of_time_array = N * (M + 1) + (N - 1) * M;
    const int N_times_M = N * M; // N is row, M is column


    const int ratio = int(size_of_time_array/5);//ratio of big step and small step
    double small_tau = big_tau/double(ratio);//small time step

    for(int i = 0 ; i < N_thread ; i++) {
        energy_integration[i] = 0.0;
    }
    //Parallel Computing start
    #pragma omp parallel num_threads(N_thread)
    {
        const int rank = omp_get_thread_num();

        double** energy_array = new double*[N];
        for(int i = 0 ; i < N ; i ++) {
          energy_array[i] = new double[M+2];
        }
        double *E_avg = new double[N_times_M];
        double* last_update = new double[N_times_M];
        for(int i = 0; i <N_times_M; i++)
        {
          E_avg[i] = 0;
          last_update[i] = 0;
        }
        trng::yarn2 r;
        trng::uniform01_dist<> u;
        r.seed(time(NULL));
        r.split(N_thread, rank);

        //Initialize energy_array
        for(int i = 0 ; i < N ; i++) {
          for(int j = 0 ; j < M + 2 ; j++) {
            if(j == 0) {
              energy_array[i][j] = TL;
            }
            else if(j == M + 1) {
              energy_array[i][j] = TR;
            }
            else {
              energy_array[i][j] = 1;
            }
          }
        }

        interaction* time_array = new interaction[size_of_time_array];
        int index = 0;
        int total_number_of_row = 2 * N - 2;
        for(int i = 0 ; i <= total_number_of_row ; i++) {
            if(i % 2 == 0) {
                for(int j = 0 ; j <= M ; j++) {
                    time_array[index].time = -log(1 - u(r))/energy_summation(energy_array, i, j, vertical); //sqrt(energy_array[index] + energy_array[index+1]);
                    time_array[index].location = index;
                    time_array[index].coordinate = make_pair(i, j);
                    time_array[index].direction = vertical;
                    time_array[index].adjacent_clocks = find_adjacent_clocks(i, j, vertical, N, M, index);
                    time_array[index].left = NULL;
                    time_array[index].right = NULL;
                    time_array[index].level = -1;
                    index++;
                }
            }
            else {
                for(int j = 0 ; j < M ; j++) {
                    time_array[index].time = -log(1 - u(r))/energy_summation(energy_array, i, j, horizontal);
                    time_array[index].location = index;
                    time_array[index].coordinate = make_pair(i, j);
                    time_array[index].direction = horizontal;
                    time_array[index].adjacent_clocks = find_adjacent_clocks(i, j, horizontal, N, M, index);
                    time_array[index].left = NULL;
                    time_array[index].right = NULL;
                    time_array[index].level = -1;
                    index++;
                }
            }
        }

        site** E_avg_profile = new site*[N];

        for(int i = 0 ; i < N ; i++) {
            E_avg_profile[i] = new site[M + 2];
        }

        //Initialization
        for(int i = 0 ; i < N ; i++) {
            for(int j = 0 ; j < M + 2 ; j++) {
                E_avg_profile[i][j].previous_time = 0;
                E_avg_profile[i][j].coordinate = make_pair(i, j);
                E_avg_profile[i][j].summation = 0;
            }
        }

        //End

        int count = 0;

        interaction** clock_time_in_step = new interaction*[ratio + 1];//each element in the array is the head of a list
        for(int i = 0; i < ratio + 1; i++)
        {
            clock_time_in_step[i] = NULL;
        }

        gettimeofday(&t1,NULL);



        for(int out_n = 0; out_n < Step; out_n++)
        {
            big_step_distribute(clock_time_in_step,time_array, size_of_time_array ,small_tau,ratio,out_n);
            for(int in_n = 0; in_n < ratio; in_n++)
            {
                if(clock_time_in_step[in_n]!= NULL)
                {
                    update(clock_time_in_step, in_n, N, M, small_tau, ratio, out_n, time_array, energy_array, u, r, &count, energy_integration[rank]);
                }
            }

            //printf("Done!\n");

            clock_time_in_step[ratio] = NULL;

        }

        gettimeofday(&t2, NULL);
        delete[] energy_array;
        delete[] E_avg;
        delete[] time_array;
        delete[] clock_time_in_step;
    }

    data<<N<<","<<M<<","<<TL<<","<<TR<<","<<current_attempt<<",";

    auto sum = 0.0;
    for(int i = 0 ; i < N_thread ; i++) {
        data<<energy_integration[i]<<", ";
        sum += energy_integration[i];
    }
    //cout<<endl;

    data<<sum/(N_thread*Step*big_tau)<<endl;
    //int size_of_energy_array = N * (M + 2); //Make the energy array in 2D
    delete[] energy_integration;

    if(current_attempt == number_of_attempts) {
      //cout<<"OK, I'm in"<<endl;
        char col = 'A' + 5 + N_thread;
        string std_div = "";
        for(int i = 1 ; i <= 6 + N_thread ; i++) {
            std_div += ",";
        }
        std_div += ("=STDEV.S(" + string(1, col) + to_string(current_line_number - number_of_attempts + 1) + ":" + string(1, col) + to_string(current_line_number) + "),");
        data<<std_div<<endl;
    }

    myfile.close();
    data.close();
    //profile.close();
}
