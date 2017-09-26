#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>

using namespace std;

#pragma once
class flowsim
{
private:
	int throughput;
	int fairness;
	int print_flows;
	int print_routers;
	int print_totals;
	int print_metrics;
	int print_metric_accuracy;
	int print_predicted;
	int tor_selection;
	int flat_selection;
	int limited_selection;
	int proposed_selection;
	int tor_metric;
	int proposed_metric;
	int bandwidth_metric;
	int eigen_metric;
	string s;
	int num_trials;
	int num_flows;
	int num_paths;
	string cap_file;
	string init_matrix;
	int trial_num;
	int aggr_max;
	int aggr_ewma;
	int aggr_mwma;
	double crash_prob;
	int debug;
	int max_eigen_iter;

	// Information about the network
	int  num_routers;
	int BW_MAX;
	int N_PEERS;
	double *bandwidth;
	int *routerType;
	double *exitProb;
	int *node_uptime;
	int *routerIndex;

	// Router Evaluation Stuff
	double rflow_bws;
	double **observed;
	double **just_observed;
	double *metrics;
	double total_metric;
	double total_bandwidth;
	double update_thresh;
	double alpha;
	int unrated;

	// Stuff involved in choosing a path
	int rand_level;
	int weight_rand_level;
	double *weights;
	double tw;
	int route_len;
	int cur_path;

	// other variables 
	int *rand_peers;
	double *rand_peers_bw;
	double *oldmetrics;
	double *row_totals;
	double *sorted_metrics;
	int *sorted_indices;

	double guard_thresh;
	//double total_metric;
	//int unrated;

	int temp;
	string line;

	// variables for the fairness
	int **flows;
	int *maxflow;
	int **routes;
	int *routesActive;
	int maxroute;
	double *sel_levels;
	double *predict;
	int *current_route;

	double *tbandwidth;
	double *flow_bws;
	double total;
	double *bandwidth_used;

	// get route variables
	int guard_selection;
	int sel;

public:
	flowsim(void);
	~flowsim(void);
	int parseCommandLine(int argc, char *argv[]);
	bool node_crash(int node_uptimte);
	int openBWCap();
	int openInitMatrix();
	int allocateArrays();
	int cleanupArrays();
	void updateMetrics();
	int simLoop();
	double drandRange(double low, double high);
	int irandRange(int low, int high);
	double median(double* myarray, int length);
	void get_route();
	void set_selection_level();
	int get_router(int pos);
	void setGuards();
	void jFairnessBW();
	void rFairnessBW();
	double combine(double a, double b);
	void cleanup();
};

