#include "flowsim.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>

using namespace std;

flowsim::flowsim(void)
{
	srand((int)time(NULL));
	// command line parameters
	throughput = 0;
	fairness = 0;
	print_flows = 0;
	print_routers = 0;
	print_totals = 0;
	print_metrics = 0;
	print_metric_accuracy = 0;
	print_predicted = 0;
	tor_selection = 0;
	flat_selection = 0;
	proposed_selection = 0;
	limited_selection = 0;
	tor_metric = 0;
	proposed_metric = 0;
	bandwidth_metric = 0;
	eigen_metric = 0;
	s = "-1";
	num_trials = 1;
	num_flows = 0;
	num_paths = 1;
	cap_file = "cap.data";
	init_matrix = "";
	trial_num = 0;
	aggr_max = 0;
	aggr_ewma = 0;
	aggr_mwma = 0;
	crash_prob = 0;
	debug = 0;
	max_eigen_iter = 0;
	
	// Information about the network
	num_routers = 0;
	BW_MAX = 10000000;
	N_PEERS = 5;

	// Stuff involved in choosing a path
	rand_level = 0;
	weight_rand_level = 0;
	
	tw = 0;
	guard_selection = 0;
	guard_thresh = 0;
	route_len = 3;
	cur_path = 0;
	crash_prob = 0;

}



flowsim::~flowsim(void)
{
}

int flowsim::parseCommandLine(int argc, char *argv[])
{
	// Parse the command line
	for(int i = 1; i < argc; i++)
	{
		if(!strcmp(argv[i],"-Dt"))
		{
			throughput = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Drf"))
		{
			fairness = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Djf"))
		{
			fairness = 2;
		}
		else if(!strcmp(argv[i],"-Pf"))
		{
			print_flows = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Pr"))
		{
			print_routers = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Pt"))
		{
			print_totals = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Pm"))
		{
			print_metrics = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Pa"))
		{
			print_metric_accuracy = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Pp"))
		{
			print_predicted = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-St"))
		{
			tor_selection = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Sf"))
		{
			flat_selection = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Sp"))
		{
			proposed_selection = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Sl"))
		{
			limited_selection = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Mt"))
		{
			tor_metric = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Mp"))
		{
			proposed_metric = 1;
			continue;
		}
		else if(!strcmp(argv[i], "-Mb"))
		{
			bandwidth_metric = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Me"))
		{
			eigen_metric = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-s"))
		{
			i++;
			if(atoi(argv[i]) == -2)
			{
				rand_level = 1;
				s = "-2";
			}
			else
			{
				// first count the number of input numbers
				temp = 0;
				while(true)
				{
					if(argv[i][0] == '-')
						break;
					temp++;
				}
				// now we need to allocate the array
				weights = new double[temp];
				s = "";
				for(int j = 0; j < temp; j++)
				{
					weights[j] = atof(argv[i+j]);
					tw += weights[j];
					s += argv[i+j];
				}
			}
			continue;
		}
		else if(!strcmp(argv[i],"-N"))
		{
			i++;
			num_trials = atoi(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-n"))
		{
			i++;
			num_flows = atoi(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-f"))
		{
			i++;
			num_paths = atoi(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-c"))
		{
			i++;
			cap_file = argv[i];
			continue;
		}
		else if(!strcmp(argv[i],"-O"))
		{
			i++;
			init_matrix = argv[i];
			continue;
		}
		else if(!strcmp(argv[i],"-Ax"))
		{
			aggr_max = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Ae"))
		{
			aggr_ewma = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-Am"))
		{
			aggr_mwma = 1;
			continue;
		}
		else if(!strcmp(argv[i],"-a"))
		{
			i++;
			alpha = atof(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-C"))
		{
			i++;
			crash_prob = atof(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-d"))
		{
			i++;
			debug = atoi(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-ei"))
		{
			i++;
			max_eigen_iter = atoi(argv[i]);
			continue;
		}
		else if(!strcmp(argv[i],"-R"))
		{
			i++;
			route_len = atoi(argv[i]);
			continue;
		}
	}

	update_thresh = tor_metric ? 2:1;
	
	if((throughput + fairness != 1) || (aggr_max + aggr_ewma + aggr_mwma != 1) || (tor_selection + flat_selection + proposed_selection != 1) || (print_metric_accuracy + print_flows + print_routers + print_totals + print_metrics + print_predicted < 1) || (tor_metric + proposed_metric + bandwidth_metric + eigen_metric != 1))
	{
		cout << "$0: calculate flow bandwidths and router utilizations\n";
		cout << "Usage:\n";
		cout << " -Dt\tdistribute bandwidth to maximize overall throughput\n";
		cout << " -Drf\tdistribute bandwidth to maximize flow fairness (Robin's)\n";
		cout << " -Djf\tdistribute bandwidth to maximize flow fairness (Josh's)\n";
		cout << " -Sf\tuse flat router selection model\n";
		cout << " -St\tuse current tor router selection model\n";
		cout << " -Sp\tuse proposed new router selection model\n";
		cout << " -Sl\tuse tor route selection with limited BW nodes\n";
		cout << " -Mt\tuse tor router metric\n";
		cout << " -Mb\tuse actual available bandwith available router metric\n";
		cout << " -Mp\tuse proposed router metric\n";
		cout << " -Me\tuse EigenRep router metric\n";
		cout << " -ei <num>\tspecify max iterations for Eigenvector power method\n";
		cout << " -Pf\tprint << the bandwidth allocated to each flow\n";
		cout << " -Pp\tprint the predicted and actual bandwidth for each flow\n";
		cout << " -Pr\tprint the bandwidth used by each router and the fractional";
		cout << " utilization\n -Pt\tprint the total bandwidths used\n";
		cout << " -Pm\tprint observation matrix\n";
		cout << " -Pa\tprint an esitmate of metric accuracy\n";
		cout << " -Ax\taggregate router measurements using the max function\n";
		cout << " -Ae\taggregate router measurements using an EWMA\n";
		cout << " -Am\taggregate router measurements using an MWMA\n";
		cout << " -C <value>\tthe crash probability for a router in each iteration\n";
		cout << " -a\tspecify alpha for EWMA (default = 0.25)\n";
		cout << " -n <num>\tset the number of flows (default: 10 per router)\n";
		cout << " -f <num>\tset the number paths per flow (default: 1)\n";
		cout << " -R <num>\tset the number of routers on each path (default: 3)\n";
		cout << " -c <file>\tset the capacity file used (default: cap.data)\n";
		cout << " -O <file>\tspecify an initial observation matrix\n";
		cout << " -N\tnumber of trials to run\n";
		cout << " -s\tspecify selection level distribution each number seperated by a space (proposed algorithm)\n";
		cout << " -d\tprint debugging information\n";
		cout << "\n";
		cout << "Note: exactly one -D, -S, -M, and -A flag and at least one -P flag";
		cout << " are required.\n";

		return 1;
	}
	return 0;
}


int flowsim::openBWCap()
{
	// Open the bandwidth cap file
	ifstream ifbw(cap_file.c_str());
	if(!ifbw.is_open())
	{
		cout << "Failed to open " << cap_file << endl;
		return 1;
	}

	temp = 0;
	while(!ifbw.eof())
	{
		getline(ifbw, line);
		temp++;
	}

	ifbw.clear();
	ifbw.seekg(0,ios::beg);

	num_routers = temp;
	total_bandwidth = 0;
	bandwidth = new double[temp];
	routerType = new int[temp];
	exitProb = new double[temp];
	routerIndex = new int[num_routers];

	
	string split;
	int index = 0;

	for(int i = 0; i < temp; i++)
	{
		getline(ifbw, line);
		index = line.find(" ");
		split = line.substr(0, index);
		line = line.substr(index + 1);
		bandwidth[i] = atof(split.c_str());
		total_bandwidth += bandwidth[i];

		index = line.find(" ");
		split = line.substr(0, index);
		line = line.substr(index + 1);
		routerType[i] = atoi(split.c_str());

		exitProb[i] = atof(line.c_str());
	}

	if(num_flows == 0)
	{
		num_flows = 1 * num_routers;
	}

	ifbw.close();
	return 0;
}

int flowsim::openInitMatrix()
{
	// Open the init matrix
	if(init_matrix != "")
	{
		ifstream imf(init_matrix);
		if(!imf.is_open())
		{
			cout << "error opening " << init_matrix << endl;
			return 1;
		}

		for(int i = 0; i < num_routers; i++)
		{
			getline(imf,line);
			if(imf.eof())
			{
				cout << "error not a valid initialization matrix not enough rows\n";
				return 1;
			}

			for(int j = 0; j < num_routers; j++)
			{
				temp = line.find(" ");
				if(temp == 0 || temp == -1)
				{
					cout << "error not a valid initialization matrix not enough columns in row " << i << endl;
					return 1;
				}
				observed[i][j] = atof(line.substr(0,temp).c_str());
				line.erase(0,temp);
			}
		}

		imf.close();
	}
	return 0;
}

int flowsim::allocateArrays()
{
	// allocate arrays
	metrics = new double[num_routers];
	observed = new double*[num_routers];
	for(int i = 0; i < num_routers; i++)
	{
		observed[i] = new double[num_routers];
	}

	just_observed = new double*[num_routers];
	for(int i = 0; i < num_routers; i++)
	{
		just_observed[i] = new double[num_routers];
	}

	bandwidth_used = new double[num_routers];
	node_uptime = new int[num_routers];
	for(int i = 0; i < num_routers; i++)
	{
		node_uptime[i] = 0;
	}

	sel_levels = new double[num_flows * num_paths];
	predict = new double[num_flows * num_paths];
	tbandwidth = new double[num_routers];
	flow_bws = new double[num_flows*num_paths];
	bandwidth_used = new double[num_routers];

	if(fairness != 0)
	{
		flows = new int*[num_routers];
		maxflow = new int[num_routers];
		for(int i = 0; i < num_routers; i++)
		{
			flows[i] = new int[num_flows*num_paths];
			maxflow[i] = 0;
		}
		
		routesActive = new int[num_flows*num_paths];
		routes = new int*[num_flows*num_paths];
		for(int i = 0; i < (num_flows*num_paths); i++)
		{
			routes[i] = new int[route_len];
		}
		maxroute = 0;

		current_route = new int[route_len];
		
		
	}

	if(proposed_metric)
	{
		rand_peers = new int[N_PEERS];
		rand_peers_bw = new double[N_PEERS];
	}

	if(eigen_metric)
	{
		oldmetrics = new double[num_routers];
		row_totals = new double[num_routers];
		sorted_indices = new int[num_routers];
		sorted_metrics = new double[num_routers];
	}
	return 1;
}

void flowsim::cleanup()
{
	// de-allocate arrays
	delete [] bandwidth;
	delete [] routerType;
	delete [] exitProb;

	delete [] metrics;

	for(int i = 0; i < num_routers; i++)
	{
		delete [] observed[i];
	}

	delete [] observed;
		
	for(int i = 0; i < num_routers; i++)
	{
		delete [] just_observed[i];
	}

	delete [] just_observed;
	delete [] bandwidth_used;

	delete [] node_uptime;

	delete [] sel_levels;
	delete [] predict;
	delete [] tbandwidth;
	delete [] flow_bws;
	delete [] bandwidth_used;

	if(fairness != 0)
	{
		for(int i = 0; i < num_routers; i++)
		{
			delete [] flows[i];
		}
		
		delete [] flows;
		delete [] maxflow;

		delete [] routesActive;
		
		for(int i = 0; i < (num_flows*num_paths); i++)
		{
			delete [] routes[i];
		}
		delete [] routes;

		delete [] current_route;
		
		
	}

	if(proposed_metric)
	{
		delete [] rand_peers;
		delete [] rand_peers_bw;
	}

	if(eigen_metric)
	{
		delete [] oldmetrics;
		delete [] row_totals;
		delete [] sorted_indices;
		delete [] sorted_metrics;
	}
	return;
}

double flowsim::drandRange(double low, double high)
{
	double random = ((double)rand())/((double)RAND_MAX);
	random = random * (high - low);
	random = random + low;
	return random;
}

int flowsim::irandRange(int low, int high)
{
	double random = ((double)rand())/((double)RAND_MAX);
	random = random * (high - low);
	random = random + low;
	return (int)random;
}

bool flowsim::node_crash(int node_uptime)
{
	double crash = ((double)rand())/((double)RAND_MAX);
	if(crash < crash_prob)
	{
		return false;
	}
	else
	{
		return true;
	}
}

void flowsim::updateMetrics()
{
	//unachievable ideal case
	if(bandwidth_metric)
	{
		for(int i = 0; i < num_routers; i++)
		{
			metrics[i] = bandwidth[i];
		}
	}

	// current case
	if(tor_metric)
	{
		for(int i = 0; i < num_routers; i++)
		{
			metrics[i] = observed[i][i];
			if(metrics[i] > BW_MAX)
				metrics[i] = BW_MAX;
		}
	}

	if(proposed_metric)
	{
		//implements the proposed metric
		for(int i = 0; i < N_PEERS; i++)
		{
			rand_peers[i] = irandRange(0, num_routers);
		}
		for(int i = 0; i < num_routers; i++)
		{
			for(int j = 0; j < N_PEERS; j++)
			{
				rand_peers_bw[j] = observed[i][rand_peers[j]];
			}
			metrics[i] = this->median(rand_peers_bw, N_PEERS);
		}
	}

	if(eigen_metric)
	{
		// get the row totals
		for(int i = 0; i < num_routers; i++)
		{
			row_totals[i] = 0;
			
			for(int j = 0; j < num_routers; j++)
			{
				row_totals[i] += observed[i][j];
			}
			if(row_totals[i] == 0)
				row_totals[i] = 1;
		}

		// initialize
		if(total_metric)
		{
			for(int i = 0; i < num_routers; i++)
			{
				oldmetrics[i] = metrics[i];
			}
		}
		else
		{
			for(int i = 0; i < num_routers; i++)
			{
				oldmetrics[i] = 1 / (double)num_routers;
			}
			if(debug)
			{
				cout << "STDERR Starting Eigenspeed from uniform vector...\n";
			}
		}

		int loopcount = 0;
		double delta = 0;

		do
		{
			for(int i = 0; i < num_routers; i++)
			{
				metrics[i] = 0;
				for(int j = 0; j < num_routers; j++)
				{
					metrics[i] += observed[j][i] * oldmetrics[j] / row_totals[j];
				}
			}

			delta = 0;
			
			for(int i = 0; i < num_routers; i++)
			{
				delta += pow(metrics[i] - oldmetrics[i],2); 
			}

			for(int i = 0; i < num_routers; i++)
			{
				oldmetrics[i] = metrics[i];
			}
			loopcount++;
		}
		while((delta > 1e-10) && (loopcount < max_eigen_iter));

		if(debug)
		{
			cout << "EigenSpeed took " << loopcount << " iterations.\n";
		}

		double swap = 0;
		total_metric = 0;
		unrated = 0;
		// now we must create the sorted_indices structure
		for(int i = 0; i < num_routers; i++)
		{
			sorted_indices[i] = i;
			sorted_metrics[i] = metrics[i];
			total_metric += metrics[i];
			if(metrics[i] == 0)
				unrated++;
		}

		for(int i = 0; i < num_routers; i++)
		{
     		for(int j = i; j < num_routers; j++)
			{
				if(sorted_metrics[j] > sorted_metrics[i])
				{
					swap = sorted_metrics[i];
					sorted_metrics[i] = sorted_metrics[j];
					sorted_metrics[j] = swap;
					sorted_indices[i] = j;
					sorted_indices[j] = i;
				}
			}
		}

		guard_thresh = metrics[sorted_indices[num_routers/2]];

		if(debug)
		{
			cout << "STDERR Unrated: " << unrated << " \n";
		}
	} // end eigen metric
}

double flowsim::median(double* myarray, int length)
{
	int zcount = 0;
	double swap;
	for(int i = 0; i < length; i++)
	{
		if(myarray[i] == 0)
			zcount++;
	}
	if(zcount > 2)
		return 0;

	// now sort the array
	for(int i = 0; i < length; i++)
	{
		for(int j = i; j < length; j++)
		{
			if(myarray[i] > myarray[j])
			{
				swap = myarray[i];
				myarray[i] = myarray[j];
				myarray[j] = swap;
			}
		}
	}

	// now return the median
	if(length % 2 == 0)
	{
		return (myarray[(int)length / 2] + myarray[((int)length / 2) + 1]) / 2.0;
	}
	else
	{
		return myarray[((int)length / 2) + 1];
	}
}

int flowsim::simLoop()
{
// Begin main loop
	while(trial_num < num_trials)
	{
		trial_num++;

		//Reset per tick stuff
		for(int i = 0; i < num_routers; i++)
		{
			bandwidth_used[i] = 0;
			for(int j = 0; j < num_routers; j++)
			{
				just_observed[i][j] = 0;
			}
		}

		for(int i = 0; i < num_flows * num_paths; i++)
		{
			flow_bws[i] = 0;
		}

		// Deal with churn
		for(int node = 0; node < num_routers; node++)
		{
			if(node_crash(node_uptime[node]))
			{
				for(int i = 0; i < num_routers; i++)
				{
					observed[node][i] = 0;
					observed[i][node] = 0;
				}
				node_uptime[node] = 0;
			}
			else
			{
				node_uptime[node]++;
			}
		}

		// then update the metrics
		this->updateMetrics();

		// allocate the flows and gather flow data
		if(throughput)
		{
			// currently not implemented
		}
		if(fairness != 0)
		{
			// reset the flows and routes indices 
			for(int i = 0; i < num_routers; i++)
				maxflow[i] = 0;
			maxroute = 0;

			for(int flow = 0; flow < (num_flows*num_paths); flow++)
			{
				this->get_route();

				for(int j = 0; j < route_len; j++)
				{
					flows[current_route[j]][maxflow[current_route[j]]] = flow;
					maxflow[current_route[j]]++;

					routes[flow][j] = current_route[j];
				}
				// Activate the flow
				routesActive[flow] = 1;
				maxroute++;

				sel_levels[flow] = sel;

				// predict the bw bottleneck of the flow
				predict[flow] = pow(2.0,30);
				for(int j = 0; j < route_len; j++)
				{
					if(metrics[routes[maxroute - 1][j]] < predict[flow])
					{
						predict[flow] = metrics[routes[maxroute - 1][j]];
					}
				}
			} // end of main flow loop

			// now run the fairness bw usage funciton
			if(fairness == 1)
			{
				this->rFairnessBW();
			}
			else if(fairness == 2)
			{
				this->jFairnessBW();
			}
		} // end fairness

		// update observation matrix
		for(int i = 0; i < num_routers; i++)
		{
			for(int j = 0; j < num_routers; j++)
			{
				if(just_observed[i][j])
				{
					observed[i][j] = this->combine(observed[i][j],just_observed[i][j]);
				}
			}
		}

		// Finally print any requested reports
		if(print_totals)
		{
			cout << total << endl;
		}

		if(print_routers)
		{
			for(int i = 0; i < num_routers; i++)
			{
				cout << bandwidth_used[i] << " " << bandwidth[i] << endl;
			}
		}

		if(print_flows)
		{
			for(int i = 0; i < num_flows * num_paths; i++)
			{
				cout << sel_levels[i] << " " << flow_bws[i] << endl;
			}
		}

		if(print_predicted)
		{
			for(int i = 0; i < num_flows * num_paths; i++)
			{
				cout << predict[i] << " " << flow_bws[i] << endl;
			}
		}

		if(print_metric_accuracy)
		{
			double t = 0;
			for(int i = 0; i < num_routers; i++)
			{
				t += ((metrics[i]/total_metric ? total_metric : 1)) - pow((bandwidth[i]/total_bandwidth),2);
			}
			cout << t << endl;
		}

	} // end main while loop

	if(print_metrics)
	{
		for(int i = 0; i < num_routers; i++)
		{
			for(int j = 0; j < num_routers; j++)
			{
				cout << observed[i][j] << " ";
			}
			cout << endl;
		}
	}

	return 0;
}

void flowsim::jFairnessBW()
{
	// make a backup of the bandwidth
	for(int i = 0; i < num_routers; i++)
	{
		tbandwidth[i] = bandwidth[i];
		bandwidth_used[i] = 0;
	}

	int end = 1;

	while(end == 1)
	{
		int min = 0;

		// get the bottleneck router
		for(int i = 0; i < num_routers; i++)
		{
			if(maxflow[i] != 0)
			{
				if((tbandwidth[i] / (double)maxflow[i]) < (tbandwidth[min] / (double)maxflow[min]))
				{
					min = i;
				}
			}
		}

		// min should now contain the index of the bottleneck router
		// increment all remaining active flows by the bottleneck amount
		double increment = tbandwidth[min] / (double)maxflow[min];
			
		for(int i = 0; i < (num_flows * num_paths); i++)
		{
			if(routesActive[i] == 1)
			{
				flow_bws[i] += increment;
			}
		}

		// remove the bottleneck bw from each router with an active flow
		for(int i = 0; i < num_routers; i++)
		{
			if(maxflow[i] > 0)
			{
				
				tbandwidth[i] -= increment * (double)maxflow[min];
				bandwidth_used[i] += increment * (double)maxflow[min];
			}
		}

		// which flows go via the bottleneck router
		// they are maxed out so remove them
		int currentFlow = 0;
		int currentRouter = 0;

		for(int i = 0; i < maxflow[min]; i++)
		{
			currentFlow = flows[min][i];
			for(int j = 0; j < route_len; j++)
			{
				currentRouter = routes[currentFlow][j];
				for (int k = 0; k < maxflow[currentRouter] - 1; k++)
				{
					if(flows[currentRouter][k] == flows[min][i])
					{
						flows[currentRouter][k] = flows[currentRouter][k+1];
						flows[currentRouter][k+1] = flows[min][i];
					}
				}
				maxflow[currentRouter]--;
			}
			// deactivate the flow
			routesActive[flows[min][i]] = 0;
			flows[min][i] = -1;
			
		}
		
		maxflow[min] = 0;
		

		// check to see if we are done
		end = 0;
		for(int i = 0; i < num_routers; i++)
		{
			if(maxflow[i] > 0)
				end = 1;
		}

	} // end while loop

	for(int i = 0; i < (num_flows * num_paths); i++)
	{
		for(int j = 0; j < route_len - 1; j++)
		{
			just_observed[routes[i][j]][routes[i][j]] += flow_bws[i];
			just_observed[routes[i][j]][routes[i][j+1]] += flow_bws[i];
			just_observed[routes[i][j+1]][routes[i][j]] += flow_bws[i];
		}
		just_observed[routes[i][route_len - 1]][routes[i][route_len - 1]] += flow_bws[i];
	}

	total = 0;
	for(int i = 0; i < num_routers; i++)
	{
		total += bandwidth_used[i];
	}

	return;
}

void flowsim::rFairnessBW()
{
	// make a backup of the bandwidth
	for(int i = 0; i < num_routers; i++)
	{
		tbandwidth[i] = bandwidth[i];
		bandwidth_used[i] = 0;
	}

	int end = 1;

	while(end == 1)
	{
		int min = 0;

		// get the bottleneck router
		for(int i = 0; i < num_routers; i++)
		{
			if(maxflow[i] != 0)
			{
				if((tbandwidth[i] / (double)maxflow[i]) < (tbandwidth[min] / (double)maxflow[min]))
				{
					min = i;
				}
			}
		}

		// min should now contain the index of the bottleneck router
		// increment all remaining active flows by the bottleneck amount
		double increment = tbandwidth[min] / (double)maxflow[min];
			
		for(int i = 0; i < (num_flows * num_paths); i++)
		{
			if(routes[i][0] != -1)
			{
				flow_bws[i] += increment;
			}
		}

		// remove the bottleneck bw from each router with an active flow
		for(int i = 0; i < num_routers; i++)
		{
			if(maxflow[i] != 0)
			{
				
				tbandwidth[i] -= increment * (double)maxflow[min];
				bandwidth_used[i] += increment * (double)maxflow[min];
			}
		}

		// which flows go via the bottleneck router
		// they are maxed out so remove them
		int currentFlow = 0;
		int currentRouter = 0;

		for(int i = 0; i < maxflow[min]; i++)
		{
			currentFlow = flows[min][i];
			for(int j = 0; j < route_len; j++)
			{
				currentRouter = routes[currentFlow][j];
				for (int k = 0; k < maxflow[currentRouter] - 1; k++)
				{
					if(flows[currentRouter][k] == flows[min][i])
					{
						flows[currentRouter][k] = flows[currentRouter][k+1];
						flows[currentRouter][k+1] = flows[min][i];
					}
				}
				maxflow[currentRouter]--;
			}
			//routes[flows[min][i]][0] = -1;
			flows[min][i] = -1;
			
		}
		
		maxflow[min] = 0;
		

		// check to see if we are done
		end = 0;
		for(int i = 0; i < num_routers; i++)
		{
			if(maxflow[i] > 0)
				end = 1;
		}

	} // end while loop

	for(int i = 0; i < (num_flows * num_paths); i++)
	{
		for(int j = 0; j < route_len - 1; j++)
		{
			just_observed[routes[i][j]][routes[i][j]] += flow_bws[i];
			just_observed[routes[i][j]][routes[i][j+1]] += flow_bws[i];
			just_observed[routes[i][j+1]][routes[i][j]] += flow_bws[i];
		}
		just_observed[routes[i][route_len - 1]][routes[i][route_len - 1]] += flow_bws[i];
	}

	total = 0;
	for(int i = 0; i < num_routers; i++)
	{
		total += bandwidth_used[i];
	}

	return;
}

// generates a new route and places it into current_route
void flowsim::get_route()
{
	guard_selection = 1;

	if(proposed_metric)
	{
		this->updateMetrics();
	}
	if(weight_rand_level || rand_level)
	{
		this->set_selection_level();
	}

	for(int i = 0; i < route_len; i++)
	{
		current_route[i] = this->get_router(i);
		guard_selection = 0;

		// we need to check to make sure we did not choose the same router
		for(int j = 0; j < i; j++)
		{
			if(current_route[j] == current_route[i])
			{
				// we chose the same router
				i--;
			}
		}

	} // end router selection loop
	return;
}

int flowsim::get_router(int pos)
{
	// no data yet just choose randomly
	if(total_metric == 0)
	{
		return this->irandRange(0, num_routers);
	}

	int router = 0;

	// flat random (no qualifications for entry position)
	if(flat_selection)
	{
		return this->irandRange(0, num_routers);
	}
	
	if(tor_selection || limited_selection)
	{
		int num_current = 0;
		double currentBW = 0;

		double exitBW = 0;
		double guardBW = 0;
		double totalBW = 0;
		double guardWeight = 0;
		double exitWeight = 0;
		double minTemp = 0;

		int router;

		if((pos != 0) && (pos != route_len - 1))
		{
			// we need the BW stats
			for(int i = 0; i < num_routers; i++)
			{
				if(routerType[i] == 2)
				{
					guardBW += metrics[i];
				}
				else if(routerType[i] == 3)
				{
					exitBW += metrics[i];
				}
				else if(routerType[i] == 4)
				{
					exitBW += metrics[i];
					guardBW += metrics[i];
				}
				totalBW += metrics[i];
			}
			// now calculate the metrics
			if(exitBW > (totalBW / 3))
			{
				exitWeight = (exitBW - (totalBW / 3)) / exitBW;
			}
			if(guardBW > (totalBW / 3))
			{
				guardWeight = (guardBW - (totalBW / 3)) / guardBW;
			}
		} // middle router stats
		
		// first we must get all the viable nodes and place them into a list
		for (int i = 0; i < num_routers; i++)
		{
			if(pos == 0)
			{
				// we are selecting a guard router
				if(routerType[i] == 2 || routerType[i] == 4)
				{
					routerIndex[num_current] = i;
					currentBW += metrics[i];
					num_current++;
				}
			}
			else if(pos == route_len - 1)
			{
				// we are selecting the exit router
				if(routerType[i] == 3 || routerType[i] == 4)
				{
					routerIndex[num_current] = i;
					currentBW += metrics[i];
					num_current++;
				}
			}
			else
			{
				// we are getting a middle router
				if(routerType[i] == 1)
				{
					routerIndex[num_current] = i;
					currentBW += metrics[i];
					num_current++;
				}
				if(routerType[i] == 2 && guardWeight != 0)
				{
					routerIndex[num_current] = i;
					currentBW += guardWeight * metrics[i];
					num_current++;
				}
				if(routerType[i] == 3 && exitWeight != 0)
				{
					routerIndex[num_current] = i;
					currentBW += exitWeight * metrics[i];
					num_current++;
				}
				if(routerType[i] == 4)
				{
					if(exitWeight != 0 && guardWeight != 0)
					{
						routerIndex[num_current] = i;
						if(guardWeight < exitWeight)
						{
							minTemp = guardWeight;
						}
						else 
						{
							minTemp = exitWeight;
						}
						routerIndex[num_current] = i;
						currentBW += minTemp * metrics[i];
						num_current++;
					}
					else if(exitWeight != 0)
					{
						routerIndex[num_current] = i;
						currentBW += exitWeight * metrics[i];
						num_current++;
					}
					else if(guardWeight != 0)
					{
						routerIndex[num_current] = i;
						currentBW += guardWeight * metrics[i];
						num_current++;
					}
				}
			}
		} // end while loop

		// now we need to sort the list
		int iSwap = 0;
		for(int i = 0; i < num_current; i++)
		{
			for(int j = i; j < num_current; j++)
			{
				if(metrics[routerIndex[j]] > metrics[routerIndex[j+1]])
				{
					iSwap = routerIndex[j];
					routerIndex[j] = routerIndex[j+1];
					routerIndex[j+1] = iSwap;
				}
			}
		}

		// we now have our list assembed we need to pick from it
		// first check to see if we need to bootstrap
		int iTemp = this->irandRange(0,num_routers);
		if(metrics[iTemp] == 0)
		{
			return iTemp;
		}
		// or do not bootstrap but use selection
		do
		{
			router = 0;
			double bw = this->drandRange(0, currentBW);
				
			while(bw > metrics[router])
			{
				bw -= metrics[router];
				router++;
			}
		}
		while((limited_selection) && (pos == route_len - 1) && (this->drandRange(0,1) > exitProb[router]));

		return router;
	}
	/*
	// current Tor 
	if(tor_selection)
	{
		if(unrated > .9*num_routers)
		{
			// this is really only the case if we are using the proposed metric 
			// because it takes a long time to bootstrap
			return this->irandRange(0, num_routers);
		}

		do
		{
			router = 0;
			double bw = this->drandRange(0, total_metric);
			while(bw > metrics[router])
			{
				bw -= metrics[router];
				router++;
			}
		}
		while (guard_selection && (metrics[router] < guard_thresh));
		return router;
	} // end current tor selection
	*/

	// proposed tor
	if(proposed_selection)
	{
		if(sel == 0)
		{
			int min = guard_selection ? num_routers/2 : 0;
			router = min + this->irandRange(0, num_routers - min);
		}
		else if(!guard_selection && (this->irandRange(0, num_routers) < unrated))
		{
			// pick randomly form the unrated routers to bootstrap them
			router = num_routers - this->irandRange(0,unrated);
		}
		else
		{
			// otherwise pick from the rated routers by selection level
			do
			{
				double dtemp = pow(2,sel * this->drandRange(0,1)) - 1;
				dtemp = dtemp / (pow(2.0, sel - 1));
				dtemp = (num_routers - unrated) * dtemp;
				router = (int)dtemp;
			}
			while(guard_selection && (router > num_routers/2));
		}
		return sorted_indices[router];
	} // end proposed selection
	return 0;
}
// we don't need this since we only are tracking the intermediate routers
/*
void flowsim::setGuards()
{
	// this function selects the guard nodes that will be used by the routers
	int num_current = 0;
	double guardBW = 0;
	// first we must get all the guard nodes and place them into a list
	for (int i = 0; i < num_routers; i++)
	{
		if(routerType[i] == 2 || routerType[i] == 4)
		{
			routerIndex[num_current] = i;
			guardBW += bandwidth[i];
			num_current++;
		}
	}

	// we now have our list of routers
	// we need to sort by their BW
	double iSwap = 0;
	for(int i = 0; i < num_current; i++)
	{
		for(int j = i; j < num_current; j++)
		{
			if(bandwidth[routerIndex[j]] > bandwidth[routerIndex[j+1]])
			{
				iSwap = routerIndex[j];
				routerIndex[j] = routerIndex[j+1];
				routerIndex[j+1] = iSwap;
			}
		}
	}

	// now we need to select the guard routers for all first routers
	

	GuardsSet = 1;
	return;
}
*/
void flowsim::set_selection_level()
{
	// all paths associated with the same flow get the same selection level
	if(cur_path < num_paths)
		return;

	cur_path = 0;
	if(weight_rand_level)
	{
		double r = this->drandRange(0,1);
		double t = weights[0];
		sel = 0;
		while(t <= r)
		{
			sel++;
			t += weights[sel];
		}
	}
	else if(rand_level)
	{
		s = this->irandRange(0,16);
	}
}

double flowsim::combine(double dold, double dnew)
{
	if(aggr_max)
	{
		if(dnew > update_thresh * dold)
		{
			return dnew;
		}
		else
		{
			return dold;
		}
	}
	if(aggr_ewma)
	{
		return (dold ? (alpha*dold + (1 - alpha)*dnew) : dnew);
	}
	if(aggr_mwma)
	{
		return (dold ? (dold + dnew) / 2 + (0.5 - alpha) * abs(dold - dnew) : dnew);
	}
	return 0;
}