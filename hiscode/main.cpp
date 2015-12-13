#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <fstream>
#include <cctype>
#include <ctime>
#include <vector>
#include <climits>

#include "mpi.h"
#include <google/dense_hash_map>
#include "pstream.h"

using namespace std;
using google::dense_hash_map;
using __gnu_cxx::hash;
using redi::pstream;

// how big should we make zero when reading output from phc
#define ZERO 0.00000000001f
#define ACTMIN 2
#define ACTMAX 10

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
  }
};

struct s_b
{
	char name[6];
	double probability;
	double imaginary;
};

int main(int argc, char *argv[])
{
	// check if arg is an input file
	if(argc != 2)
	{
		fprintf(stderr, "Usage: nplayer_se [input.game]\n");
		return 0;
	}
	
	// mpi vars
	int n_id, n_sz; // this processor id, number of processors
	MPI_Init(&argc, &argv); // start your engines
	MPI_Comm_rank(MPI_COMM_WORLD, &n_id); // get this processor id
	MPI_Comm_size(MPI_COMM_WORLD, &n_sz); // number of processors
	MPI_Status st; // status var for recieve
	
	// other vars
	int i, j, k, l; // iterators, temp vars
	int* z; // in case we need a pointer
	fstream fp; // file pointer for data input
	int p_sz; // number of players
	int* p_act; // number of actions per player
	int** ch_a; // for help enumerating over actions
	int* p_k; // combination size per player
	int p_act_min; // minimum number of actions out of all players
	
	dense_hash_map<const char*, int*, hash<const char*>, eqstr> p; // payoff hash map
	dense_hash_map<const char*, int*, hash<const char*>, eqstr>::iterator it; // payoff hash map
	char* k_b;
	
	//map<string, int*> p;
	//map<string, int*>::iterator it;
	
	string kval (" "); // for access to hash map
	string input_buffer, temp_buffer, o_b; // input buffer for line input, temp for help with substring, output buffer
	char t_b[100]; // yet another temporary buffer
	clock_t clk, clk_part; // timer
	int done, next; // doneness for permutation iterator
	int* ch_j; // for chase algorithm
	int* ch_r; // ..
	int** ch_w; // ..
	unsigned long n_comb; // running tally of combinations
	int n_comb_s;
	
	pstream phc; // for running phc
	string phc_s, eq_i, eq_t, key_b; // ..
	int phc_unk, phc_eq; // number of unknowns
	
	int** eq_m; // for checking combinations for NE
	int* eq_c; // ..
	int eq_d; // .. even more
	int eq_0;
	string eq_buf[10]; // equation buffer
	
	s_b* sol_b; // solution buffer
	int sol_f; // solution found
	string r_b; // result buffer
	
	// timer start
	clk = clock();
	double t0 = MPI_Wtime();
	// open input
	fp.open(argv[1], fstream::in);
	
	if(!fp.is_open())
	{
		fprintf(stderr, "Cannot open input file.\n");
		return 0;
	}
	
	// get number of players, create one dimension
	getline(fp, input_buffer);
	p_sz = atoi(input_buffer.c_str());
	input_buffer.clear();
	
	// create space for arrays for number of actions for each player and enumerating through each
	p_act = (int*)malloc(p_sz * sizeof(int));
	ch_a = (int**)malloc(p_sz * sizeof(int*));
	p_k = (int*)malloc(p_sz * sizeof(int));
	ch_j = (int*)malloc(p_sz * sizeof(int));
	ch_r = (int*)malloc(p_sz * sizeof(int));
	ch_w = (int**)malloc(p_sz * sizeof(int*));
	eq_c = (int*)malloc(p_sz * sizeof(int*));
	
	// get number of actions for each player -- lets hope there's no more than 30 or so players
	getline(fp, input_buffer);
	
	// iterate through it
	j = 0;
	k = 0;
	p_act_min = INT_MAX;
	for(i = 0; k < p_sz; i++)
	{
		// we have edge of number, throw it is in to array
		if(!isdigit(input_buffer[i]) && input_buffer[i] != '-')
		{
			// number of actions per player
			temp_buffer = input_buffer.substr(j, i - j);
			p_act[k] = atoi(temp_buffer.c_str());
			
			// check if a player only has less than 2 actions
			if(p_act[k] < ACTMIN)
			{
				cout << "Player " << k << " has < 2 pure strategies.  Exiting." << endl;
				return 0;
			}
			
			// check if this is a minimum
			if(p_act_min > p_act[k]) p_act_min = p_act[k];
			
			// make enumerator, helper for chase algo
			ch_a[k] = (int*)malloc(p_act[k] * sizeof(int*));
			ch_w[k] = (int*)malloc(p_act[k] * sizeof(int*));
			
			// next
			k++;
			j = i;
		}
	}
	
	// init hash map
	p.set_empty_key(NULL);
	
	// start loading crap in to hash map
	while(!fp.eof())
	{
		// clear hash key and other strings
		kval.clear();
		input_buffer.clear();
		temp_buffer.clear();
		
		// get next line
		getline(fp, input_buffer);
		
		// make sure we have actual input
		if(input_buffer.length() > 0)
		{
			// iterate through it, getting whatever is between [...]
			j = 0;
			for(i = 1; i < input_buffer.length(); i++)
			{
				// find the first ]
				if(input_buffer[i] == ']')
					break;
				j++;
			}
			
			// make kval -- index, comma, repeat, ends in \0 because of .c_str()
			temp_buffer = input_buffer.substr(1, i - 1);
			for(j = 0; j < temp_buffer.length(); j++)
			{
				if(temp_buffer[j] != 32)
				{
					kval += temp_buffer[j];
					kval += ",";
				}
			}
			
			// make new space for it, +1 for null termination
			k_b = (char*)malloc(kval.length() + 1);
			
			// copy in to k_b
			strcpy(k_b, kval.c_str());
			
			// check to see if key has been used
			it = p.find(k_b);
			if(it != p.end())
			{
				// found key, exist
				cout << "Found duplicate key: " << k_b << ".  Exiting." << endl;
				return 0;
			}
			
			// use kval has a key to input payoffs, create space
			p[k_b] = (int*)malloc(p_sz * sizeof(int));
			
			// skip characters until we have found the next digit or - sign
			while(!isdigit(input_buffer[++i]));
			if(input_buffer[i - 1] == '-') i--;
		
			// iterate through the next p_sz numbers in this line
			j = i;
			k = 0;
			for(; k < p_sz; i++)
			{	
				// we have edge of number, throw it is in to array
				if(!isdigit(input_buffer[i]) && input_buffer[i] != '-')
				{
					temp_buffer = input_buffer.substr(j, i - j);
					p[k_b][k] = atoi(temp_buffer.c_str());
					k++;
					j = i;
				}
			
				// clear the temp buffer
				temp_buffer.clear();
			}
		}
	}
	
	// close file pointer
	fp.close();
	
	/*
	// dump values
	for(it = p.begin(); it != p.end(); it++)
	{
		cout << (*it).first << " -> " << (*it).second << endl;
	}
	*/
	
	// data is now loaded -- see how long it took
	
	cout << "Node " << n_id << ": Loaded data in " << MPI_Wtime()-t0 << " seconds\n";
	
	// timer
	clk = clock();
	double t1 = MPI_Wtime();
	// reset permutation for k = 2 for all players
	k = ACTMIN;
	for(i = 0; i < p_sz; i++)
	{
		for(j = 0; j < p_act[i]; j++)
		{
			// a[j] <- 0 for 0 <= j < s, a[j] <- 1 for s <= j < n
			j < (p_act[i] - k) ? ch_a[i][j] = 0 : ch_a[i][j] = 1;
			
			// w[j] <- 1 for 0 <= j < n
			ch_w[i][j] = 1;
		}
		
		// if s > 0, set r <- s, otherwise set r <- t
		p_act[i] - k > 0 ? ch_r[i] = p_act[i] - k : ch_r[i] = k;
		p_k[i] = k;
	}
	
	// now permutate
	done = false;
	n_comb = 0;
	n_comb_s = 0;
	clk_part = clock();
	do
	{
		// increment combinations
		n_comb++;
		
		if(n_comb % n_sz == n_id)
		{
			n_comb_s++;
			if(n_comb_s % 1000 == 0)
			{
				cout << "Node " << n_id << " chewed thorugh 1000 combinations in " << ((double)(clock() - clk_part)) / CLOCKS_PER_SEC << " seconds.\n";
				clk_part = clock();
			}
			
			/*
			sleep(2);
			// we have a combination, print it
			cout << "Node " << n_id << ": {";
			l = 0;

			for(i = 0; i < p_sz; i++)
			{
				cout << "(";
				o_b.clear();
				l = 0;

				for(j = 0; j < p_act[i]; j++)
				{
					if(ch_a[i][j] == 1)
					{
						sprintf(t_b, "%d", j + 1);
						o_b += t_b;
						o_b += ", ";
						l++;
					}
				}

				cout << o_b.substr(0, (p_k[i] * 3) - 2); // "1, 2, 3, 4, 5, "

				if(l > p_k[i]){ cout << ", ERROR " << l << " < " << p_k[i] << endl; return 0; }

				i + 1 == p_sz ? cout << ")" : cout << "), ";
			}
			cout << "}\n";
			*/

		
			// debugging nonsense for chase
			/*
			for(i = 0; i < p_sz; i++)
			{
				cout << "player " << i << "\n";
				cout << "ch_j = " << ch_j[i] << "\n";
				cout << "ch_r = " << ch_r[i] << "\n";
				cout << "ch_w = {";
			
				for(j = 0; j < p_act[i]; j++)
					j == p_act[i] - 1 ? cout << ch_w[i][j] << "}\n" : cout << ch_w[i][j] << ", ";
				
				cout << "a = {";
			
				for(j = 0; j < p_act[i]; j++)
					j == p_act[i] - 1 ? cout << ch_a[i][j] << "}\n" : cout << ch_a[i][j] << ", ";
				
				cout << "\n\n";
			}
			*/
		

			// remake the map
			eq_m = (int**)malloc(p_sz * sizeof(int*));
		
			for(i = 0; i < p_sz; i++)
			{
				eq_m[i] = (int*)malloc(p_k[i] * sizeof(int*));
				eq_c[i] = 0;
				k = 0;
			
				for(j = 0; j < p_k[i]; j++)
				{
					// find the next 1 in a
					while(ch_a[i][k] == 0) k++;
				
					// throw it in to the map
					eq_m[i][j] = k + 1;
				
					// go to next
					k++;
				}
			}
		
			// for debugging the map
			/*
			cout << "MAP" << endl;
			for(i = 0; i < p_sz; i++)
			{
				cout << "player " << i;
				for(j = 0; j < p_k[i]; j++)
				{
					cout << " " << j << " -> " << eq_m[i][j] << ";";
				}
				cout << endl;
			}
			cout << endl;
			*/
		
			// count number of unknowns
			phc_unk = 0;
			for(i = 0; i < p_sz; i++) phc_unk += p_k[i];
		
			// clear the system of equations
			phc_s.clear();
		
			// for each player check nash condition
			phc_eq = 0;
			for(i = 0; i < p_sz; i++)
			{
				/*
				cout << "now on player " << i << endl;
				*/
			
				// for each action in p_k[i]
				for(k = 0; k < p_k[i]; k++)
				{
					// clear the equation for player i action k
					eq_i.clear();

					// while there are combinations left
					eq_d = false;
					do
					{
						/*
						// print out this combination
						for(j = 0; j < p_sz; j++) cout << eq_c[j] << ", ";
						cout << endl;
						*/
					
						// clear key, reset var for checking if term is 0
						key_b.clear();
						eq_0 = false;
					
						// build key
						for(j = 0; j < p_sz; j++)
						{
							if(i != j)
								sprintf(t_b, "%d,", eq_m[j][eq_c[j]]);
							else
								sprintf(t_b, "%d,", eq_m[j][k]);
						
							key_b += t_b;
						}

						// clear this term
						eq_t.clear();

						// for each player
						for(j = 0; j < p_sz; j++)
						{
							if(eq_t.length()) eq_t += "* ";

							// player is the player we're checking NE
							if(j == i)
							{
								if(p[key_b.c_str()][j] == 0) eq_0 = true;
								sprintf(t_b, "%d ", p[key_b.c_str()][j]);
								eq_t += t_b;
							}
							// not the player we're checking NE
							else
							{
								// add a "* x## " based on player and action
								sprintf(t_b, "x%d%d ", j + 1, eq_m[j][eq_c[j]]);
								eq_t += t_b;
							}

							/*
							// output the term
							cout << eq_t << endl;
							*/
						}

						// finally add term add partial equation to player's equation in the array
						if(eq_buf[k].length() == 0) eq_t = eq_t;
						else
						{
							if(eq_t[0] == '-') eq_t = "- " + eq_t.substr(1);
							else eq_t = "+ " + eq_t;
						}

						if(eq_0 == false)
							eq_buf[k] += eq_t;

						// make next combination
						for(j = 0; j < p_sz; j++)
						{
							// skip player i
							if(i == j)
							{
								if(j + 1 == p_sz) eq_d = true;
							
								continue;
							}
						
							// we're not done incrementing this player
							if(eq_c[j] < p_k[j] - 1)
							{
								// go to next
								eq_c[j]++;

								// stop
								break;
							}
							// done with this player
							else
							{
								// reset this player
								eq_c[j] = 0;

								// test if we're done
								if(j + 1 == p_sz)
								{ 
									eq_d = true; 
									break;
								}
							}
						}
					} while(!eq_d);
				
					// reset combinations
					for(j = 0; j < p_sz; j++) eq_c[j] = 0;
				}
			
				// make combinations for phc solver .. for equations 2 (read: 1) through k_i
				for(k = 1; k < p_k[i]; k++)
				{
					// add to system: eq 1 - eq k = 0
					phc_s += eq_buf[0] + "- (" + eq_buf[k] + ");\n";
					
					phc_eq++;
				}
			
				/*
				cout << "\n\n\n\n\n";
				// what do have so far
				for(j = 0; j < 10; j++)
				{
					cout << j << " -> "<< eq_buf[j] << endl << endl;
				}
				*/

			
				/*
				// add player's equation to the system
				sprintf(t_b, " - u%d;\n", i);
				phc_s += eq_i + t_b;
				*/
			
				/*
				// output equation so far
				cout << eq_b << endl;
				*/
				// clear out everything
				for(k = 0; k < 10; k++) eq_buf[k].clear();
			}
		
			// now add the probability equations
			phc_unk = 0;
			for(i = 0; i < p_sz; i++)
			{
				// add another equation, clear some stuff
				eq_i.clear();
			
				// add each mixed strategy proability
				for(j = 0; j < p_k[i]; j++)
				{
					if(eq_i.length()) eq_i += "+ ";
					sprintf(t_b, "x%d%d ", i + 1, eq_m[i][j]);
					eq_i += t_b;
					
					phc_unk++;
				}
			
				// add equals 1
				phc_s += eq_i + "- 1;\n";
				phc_eq++;
			}
		
			// output what we have
			//cout << phc_s << endl;
			//cout << phc_eq << " equations\n";
			//cout << phc_unk << " unknowns\n";
		
			// clear the map
			for(i = 0; i < p_sz; i++) free(eq_m[i]);
		
			// call PHCpack and make sure it's open
			phc.open("./phc -b");
			if(!phc.is_open())
			{
				cout << "Node " << n_id << ": Could not open the phc pack executable.\n";
			}
		
			// start inputting
			phc << "n" << endl; // is system on input file? no
			phc << phc_eq << endl; // number of polynomials
			phc << phc_unk << endl; // number of unknowns
			
			// while there are still characters to output
			while(phc_s.length() > 0)
			{
				// get the first 100 characters
				temp_buffer = phc_s.substr(0, 100);
				phc_s.erase(0, 100);
				phc << temp_buffer << endl;
			}
			
			phc << "n" << endl; // no output file
		
			// output execution
			k = 0;
			
			while(phc >> o_b)
			{
				// we might be getting somewhere
				if(o_b == "solution")
				{	
					phc >> o_b;
				
					if(o_b == "for")
					{
						//cout << "Found a solution\n";
						k++;
						
						// ok we have a solution, skip next 2 lines
						phc >> o_b; // t
						phc >> o_b; // :
						phc >> o_b; // another :?
					
						// create space for solutions
						sol_b = (s_b*)malloc(phc_unk * sizeof(s_b));
					
						// iterate through solutions, which come in the format:
						// x(player)(pure strategy)
						// :
						// float - probability in mixed strategy
						// float - imaginary, which we hope is close to 0
						for(i = 0; i < phc_unk; i++)
						{
							//cout << "var name: " << o_b << endl;
						
							// put variable name in
							strncpy(sol_b[i].name, o_b.c_str(), 6);
							phc >> o_b;
						
							// skip :
							//cout << "a : " << o_b << endl;
							phc >> o_b;
						
							//cout << "prob: " << o_b << endl;
							// probability
							sol_b[i].probability = atof(o_b.c_str());
							phc >> o_b;
						
							//cout << "imaginary: " << o_b << endl << endl;
							// imaginary poriton
							sol_b[i].imaginary = atof(o_b.c_str());
							phc >> o_b;
						}
					
						// check solution for validity
						sol_f = true;
						for(i = 0; i < phc_unk; i++)
						{
							// make sure probability is positive
							if(sol_b[i].probability <= ZERO || sol_b[i].probability > 1) sol_f = false;
						
							// make sure imaginary is as close to 0 as possible
							if(sol_b[i].imaginary > ZERO) sol_f = false;
						}
					
						// print out results
						if(sol_f)
						{
							// clear buffer
							r_b.clear();
						
							// found solution, print combination
							sprintf(t_b, "Node %d: found solution for {", n_id);
							r_b += t_b;
							l = 0;

							for(i = 0; i < p_sz; i++)
							{
								r_b += "(";
								o_b.clear();
								l = 0;

								for(j = 0; j < p_act[i]; j++)
								{
									if(ch_a[i][j] == 1)
									{
										sprintf(t_b, "%d", j + 1);
										o_b += t_b;
										o_b += ", ";
										l++;
									}
								}

								r_b += o_b.substr(0, (p_k[i] * 3) - 2); // "1, 2, 3, 4, 5, "

								if(l > p_k[i]){ cout << ", ERROR " << l << " < " << p_k[i] << endl; return 0; }

								i + 1 == p_sz ? r_b += ")" : r_b += "), ";
							}
							r_b += "}\n";
						
							// print probabilities						
							for(i = 0; i < phc_unk; i++)
							{
								sprintf(t_b, "%s -> %f", sol_b[i].name, sol_b[i].probability);
								r_b += t_b;
								if(i + 1 != phc_unk) r_b += ", "; 
							}
							r_b += "\n";
						
							// output
							cout << r_b;
						}
					
						// free memory
						free(sol_b);
					}
				}
			}
			
			// close the process
			phc.clear();
			phc.close();
		}
		
		// generate next combinations from left to right, rightmost player being last
		next = true;
		for(i = 0; next && i < p_sz;)
		{	
			next = false;
		
			// find j and branch -- set j <- r
			ch_j[i] = ch_r[i];

			// repeat until w[j] = 1
			while(ch_w[i][ch_j[i]] == 0 && p_k[i] < p_act[i])
			{
				// set w[j] <- 1, j <- j + 1
				ch_w[i][ch_j[i]] = 1;
				ch_j[i] = ch_j[i] + 1;

				// terminate if j = n
				if(ch_j[i] == p_act[i])
				{
					// this player has no more combinations of size p_k[i] anymore
					next = true;
				
					// terminate
					break;
				}
			}
		
			// no more to make for this player of this size
			if(next || p_k[i] == p_act[i])
			{	
				// check if we can make combinations bigger
				if(p_k[i] < p_act[i])
				{
					// we can, so lets do
					p_k[i] = p_k[i] + 1;
				}
				else
				{
					// player reached its limit, reset size of combination
					p_k[i] = ACTMIN;
				}
			
				// reset this player for next size of combination
				for(j = 0; j < p_act[i]; j++)
				{
					// make it 0^s1^t
					// a[j] <- 0 for 0 <= j < s, a[j] <- 1 for s <= j < n
					j < (p_act[i] - p_k[i]) ? ch_a[i][j] = 0 : ch_a[i][j] = 1;

					// w[j] <- 1 for 0 <= j < n
					ch_w[i][j] = 1;
				}
			
				// if s > 0, set r <- s, otherwise set r <- t
				p_act[i] - p_k[i] > 0 ? ch_r[i] = p_act[i] - p_k[i] : ch_r[i] = p_k[i];
			
				// did we reset this player's size? if so go to next player
				if(p_k[i] == ACTMIN)
				{
					i++;
					next = true;
				}
				else next = false;
			
				// here we go again
				continue;
			}
		
			// did not terminate, set w[j] <- 0
			ch_w[i][ch_j[i]] = 0;

			// 4-way branch -- (c4) j is odd and a[j] not equal to 0
			if(ch_j[i] % 2 == 1 && ch_a[i][ch_j[i]] != 0)
			{
				l = 4;
			}
			// (c5) j is even and a[j] not equal to 0
			else if(ch_j[i] % 2 == 0 && ch_a[i][ch_j[i]] != 0)
			{
				l = 5;
			}
			// (c6) j is even and a[j] is 0
			else if(ch_j[i] % 2 == 0 && ch_a[i][ch_j[i]] == 0)
			{
				l = 6;
			}
			// (c7) j is odd and a[j] is 0
			else //if(j % 2 == 0 && a[j] == 0) -- shouldn't have to check for this, only valid option left
			{
				l = 7;
			}
		
			switch(l)
			{
				case 7: // c7 - move left 2

				// if a[j - 1] not equal 0, go to c6 -- if a[j - 1] is 0, do not go to c6
				if(ch_a[i][ch_j[i] - 1] == 0)
				{
					// set a[j] <- 1, a[j - 2] <- 0
					ch_a[i][ch_j[i]] = 1;
					ch_a[i][ch_j[i] - 2] = 0;

					// if r = j - 2
					if(ch_r[i] == ch_j[i] - 2)
						// set r <- j
						ch_r[i] = ch_j[i];
					// otherwise if r = j - 1
					else if(ch_r[i] == ch_j[i] - 1)
						// set r <- j - 2
						ch_r[i] = ch_j[i] - 2;

					// do not go to c6
					break;
				}

				case 6: // c6 - move left 1

				// set a[j] <- 1, a[j - 1] <- 0
				ch_a[i][ch_j[i]] = 1;
				ch_a[i][ch_j[i] - 1] = 0;

				// if r = j > 1
				if(ch_r[i] == ch_j[i] && ch_j[i] > 1)
					// set r <- j - 1
					ch_r[i] = ch_j[i] - 1;
				// otherwise if r = j - 1
				else if(ch_r[i] == ch_j[i] - 1)
					// set r <- j
					ch_r[i] = ch_j[i];

				break;

				case 5: // c5 - move right 2

				// if a[j - 2] not 0, go to 4 -- or if a[j - 2] is 0, do not go to 4
				if(ch_a[i][ch_j[i] - 2] == 0)
				{
					// otherwise set a[j - 2] <- 1, a[j] <- 0
					ch_a[i][ch_j[i] - 2] = 1;
					ch_a[i][ch_j[i]] = 0;

					// if r = j
					if(ch_r[i] == ch_j[i])
					{
						// set r <- max(j - 2, 1)
						if(ch_j[i] - 2 > 1)
							ch_r[i] = ch_j[i] - 2;
						else
							ch_r[i] = 1;
					}
					// otherwise if r = j - 2
					else if(ch_r[i] == ch_j[i] - 2)
					{
						// set r <- j - 1
						ch_r[i] = ch_j[i] - 1;
					}

					// do not go to c4
					break;
				}

				case 4: // c4 - move right 1

				// move right one, set a[j - 1] <- 1, a[j] <- 0
				ch_a[i][ch_j[i] - 1] = 1;
				ch_a[i][ch_j[i]] = 0;

				// if r = j > 1
				if(ch_r[i] == ch_j[i] && ch_j[i] > 1)
					// set r <- j - 1
					ch_r[i] = ch_j[i] - 1;
				// otherwise if r = j - 1
				else if(ch_r[i] == ch_j[i] - 1)
					// set r <- j
					ch_r[i] = ch_j[i];

				break;
			}
		
			// end chase algo
		}
	
		// rightmost player has turned over, we're done
		if(next && i == p_sz) done = true;
	
		// check if done
	} while(!done);
	double t2 = MPI_Wtime();
	cout << "Node " << n_id << ": Ran " << n_comb_s << " combinations in " << t2-t1 << " seconds\n";
	
	// clean up MPI
	if(n_sz > 1) MPI_Finalize();

	// end program
	return 1;
}
