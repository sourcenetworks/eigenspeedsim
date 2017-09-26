#!/usr/bin/perl -w
#
#use Math::MatrixSparse;
#use Math::GLPK qw(:constants);
use strict;
$| = 1; # non-buffered I/O

#command line params
my $fairness = 0;
my $throughput = 0;
my $print_flows = 0;
my $print_routers = 0;
my $print_totals = 0;
my $print_metrics = 0;
my $print_metric_accuracy = 0;
my $print_predicted = 0;
my $flat_selection = 0;
my $tor_selection = 0;
my $proposed_selection = 0;
my $tor_metric = 0;
my $proposed_metric = 0;
my $bandwidth_metric = 0;
my $eigen_metric = 0;
my $s = -1;
my $cap_file = "cap.data";
my $init_matrix = "";
my $num_flows;
my $num_trials = 1;
my $trial_num = 0;
my $aggr_max = 0;
my $aggr_ewma = 0;
my $aggr_mwma = 0;
my $crash_prob = 0;
my $debug = 0;
my $max_eigen_iter = 0;
my $num_paths = 1;

# Information about the network
my $num_routers; #determined by length of capacity file
my $BW_MAX = 10000000;
my $N_PEERS = 5; #how many network status documetns do we collect?
my @bandwidth; #real available bandwidth
my @node_uptime;

# router evaluation stuff
my ($total, @bandwidth_used, @flow_bws, $rflow_bws);
my @observed; #bandwidth tracked by peers and self
my @just_observed; #ditto, but just this iteration
my @predict; # what bandwidth do I think I'll get from this flow?
my @metrics;
my @sorted_indices;
my $total_metric = 0;
my $total_bandwidth = 0;
my $update_thresh; #tor doesn't push an update until a factor-of-2 change
my @sel_levels; #track which selection level a flow is at
my $alpha = 0.25; # for EWMA
my $unrated;


# stuff involved in choosing a path
my $rand_level = 0;
my $weight_rand_level = 0;
my @weights;
my $guard_selection = 0; # are we currently choosing a guard node?
my $guard_thresh;
my $route_len = 3;
my $cur_path = 0;

while($_ = shift)
{
	$throughput = 1 if $_ eq "-Dt";
	$fairness = 1 if $_ eq "-Df";
	$print_flows = 1 if $_ eq "-Pf";
	$print_routers = 1 if $_ eq "-Pr";
	$print_totals = 1 if $_ eq "-Pt";
	$print_metrics = 1 if $_ eq "-Pm";
	$print_metric_accuracy = 1 if $_ eq "-Pa";
	$print_predicted = 1 if $_ eq "-Pp";
	$tor_selection = 1 if $_ eq "-St";
	$flat_selection = 1 if $_ eq "-Sf";
	$proposed_selection = 1 if $_ eq "-Sp";
	$tor_metric = 1 if $_ eq "-Mt";
	$proposed_metric = 1 if $_ eq "-Mp";;
	$bandwidth_metric = 1 if $_ eq "-Mb";
	$eigen_metric = 1 if $_ eq "-Me";
	$s = shift if $_ eq "-s";
	$num_trials = shift if $_ eq "-N";
	$num_flows = shift if $_ eq "-n";
	$num_paths = shift if $_ eq "-f";
	$cap_file = shift if $_ eq "-c";
	$init_matrix = shift if $_ eq "-O";
	$aggr_max = 1 if $_ eq "-Ax";
	$aggr_ewma = 1 if $_ eq "-Ae";
	$aggr_mwma = 1 if $_ eq "-Am";
	$alpha = shift if $_ eq "-a";
	$crash_prob = shift if $_ eq "-C";
	$debug = 1 if $_ eq "-d";
	$max_eigen_iter = shift if $_ eq "-ei";
	$route_len = shift if $_ eq "-R";
}
if($s =~ /,/)
{
	$weight_rand_level = 1;
	@weights = split ',', $s;
	my $tw = 0;
	$tw += $_ for @weights;
	$weights[$_] /= $tw for (0..$#weights);
}
elsif($s == -2)
{
	$rand_level = 1;
}
$update_thresh = $tor_metric ? 2 : 1;

if(($throughput + $fairness != 1) or
   ($aggr_max + $aggr_ewma + $aggr_mwma != 1) or
   ($tor_selection + $flat_selection + $proposed_selection != 1) or
   ($print_metric_accuracy + $print_flows + $print_routers +
    $print_totals + $print_metrics + $print_predicted < 1) or
   ($tor_metric + $proposed_metric + $bandwidth_metric + $eigen_metric != 1))
{
	print "$0: calculate flow bandwidths and router utilizations\n";
	print "Usage:\n";
	print " -Dt\tdistribute bandwidth to maximize overall throughput\n";
	print " -Df\tdistribute bandwidth to maximize flow fairness\n";
	print " -Sf\tuse flat router selection model\n";
	print " -St\tuse current tor router selection model\n";
	print " -Sp\tuse proposed new router selection model\n";
	print " -Mt\tuse tor router metric\n";
	print " -Mb\tuse actual available bandwith available router metric\n";
	print " -Mp\tuse proposed router metric\n";
	print " -Me\tuse EigenRep router metric\n";
	print " -ei <num>\tspecify max iterations for Eigenvector power method\n";
	print " -Pf\tprint the bandwidth allocated to each flow\n";
	print " -Pp\tprint the predicted and actual bandwidth for each flow\n";
	print " -Pr\tprint the bandwidth used by each router and the fractional";
	print " utilization\n -Pt\tprint the total bandwidths used\n";
	print " -Pm\tprint observation matrix\n";
	print " -Pa\tprint an esitmate of metric accuracy\n";
	print " -Ax\taggregate router measurements using the max function\n";
	print " -Ae\taggregate router measurements using an EWMA\n";
	print " -Am\taggregate router measurements using an MWMA\n";
	print " -C <value>\tthe crash probability for a router in each iteration\n";
	print " -a\tspecify alpha for EWMA (default = 0.25)\n";
	print " -n <num>\tset the number of flows (default: 10 per router)\n";
	print " -f <num>\tset the number paths per flow (default: 1)\n";
	print " -R <num>\tset the number of routers on each path (default: 3)\n";
	print " -c <file>\tset the capacity file used (default: cap.data)\n";
	print " -O <file>\tspecify an initial observation matrix\n";
	print " -N\tnumber of trials to run\n";
	print " -s\tspecify selection level distribution (proposed algorithm)\n";
	print " -d\tprint debugging information\n";
	print "\n";
	print "Note: exactly one -D, -S, -M, and -A flag and at least one -P flag";
	print " are required.\n";

	exit 1;
}

open BW, $cap_file or die "Could not open $cap_file: $!\n";
@bandwidth = <BW>;
close BW;
chomp @bandwidth;

# Do initialization stuff
$num_routers = @bandwidth;
$total_bandwidth += $bandwidth[$_] for (0..$num_routers-1);
$num_flows = 10 * $num_routers unless $num_flows;

open(INIT_MATRIX, $init_matrix) if($init_matrix);
for my $i (0 .. $num_routers-1)
{
	if($init_matrix)
	{
		$_ = <INIT_MATRIX>;
		$observed[$i] = [split];
		die "Error: $init_matrix isn't a valid observation matrix.\n"
			unless @{$observed[$i]} == $num_routers;
		#print STDERR $observed[$i][$_] for (0 .. $num_routers-1);
		#print "\n";
	}
	else
	{
		$observed[$i][$_] = 0 for (0 .. $num_routers-1);
	}
}
$_ = <INIT_MATRIX> if $init_matrix;
die "Error: $init_matrix isn't a valid observation matrix" if defined $_;
close(INIT_MATRIX) if ($init_matrix);

while($trial_num++ < $num_trials)
{
	#Reset per-tick stuff
	for my $i (0 .. $num_routers-1)
	{
		$bandwidth_used[$i] = 0;
		$just_observed[$i][$_] = 0 for (0 .. $num_routers-1);
	}

	#deal with churn (done here in case observed is pre-initialized)
	for my $node (0..$num_routers-1)
	{
		if(&node_crash($node_uptime[$node]))
		{
			$observed[$node][$_]=$observed[$_][$node]=0 for (0..$num_routers-1);
			$node_uptime[$node] = 0;
		}
		else
		{
			$node_uptime[$node]++;
		}
	}
	# and then update the metrics
	&update_metrics();

	# Allocate flows and gather flow info
	if($throughput)
	{
		my $objective;
		my $flows;
		push @$flows, [ &get_route_col() ] for (1 .. $num_flows * $num_paths);
		for my $i (0 .. $num_routers - 1)
		{
			$flows->[($num_flows * $num_paths)+$i] = [];
			push @{$flows->[($num_flows*$num_paths)+$i]}, 0 for (0..$#bandwidth);
			$flows->[($num_flows*$num_paths)+$i][$i] = 1;
		}
		push @$objective, 1 for (1..$num_flows * $num_paths);
		push @$objective, 0 for (1..$num_routers);

		($total, $rflow_bws) = &throughput($flows, $objective, @bandwidth);
	}
	if($fairness)
	{
		my @flows;
		my @routes;
		$flows[$num_routers - 1] = undef;
		foreach my $flow (0 .. ($num_flows*$num_paths)-1)
		{
			foreach my $host (&get_route)
			{
				push @{$flows[$host]}, $flow;
				push @{$routes[$flow]}, $host; 
			}
			$sel_levels[$flow] = $s;
			$predict[$flow] = &min(map {$metrics[$_]} @{$routes[$flow]});
		}
		($total, $rflow_bws) = &fairness(\@routes, \@flows, \@bandwidth);
	}
	@flow_bws = @$rflow_bws;

	#update observation matrix
	for my $i (0..$num_routers-1)
	{
		for my $j (0..$num_routers-1)
		{
			$observed[$i][$j] =
					&combine($observed[$i][$j],$just_observed[$i][$j])
					if $just_observed[$i][$j];
		}
	}

	# Finally, print any requested reports
	print "$total\n" if $print_totals;

	if($print_routers)
	{
		print "$bandwidth_used[$_] $bandwidth[$_]\n" for (0..$num_routers-1);
	}

	if($print_flows)
	{
		#XXX collapse multiflow stuff here once we know it works
		for my $flow (0..$num_flows*$num_paths-1)
		{
			print "$sel_levels[$flow] $flow_bws[$flow]\n";
		}
	}
	if($print_predicted)
	{
		for my $flow (0..$num_flows*$num_paths-1)
		{
			print "$predict[$flow] $flow_bws[$flow]\n";
		}
	}

	if($print_metric_accuracy)
	{
		#print the norm of @bandwidth/$total_bandwidth-@metrics/$total_metrics
		my $t = 0;
		$t += (($metrics[$_]/($total_metric ? $total_metric : 1)) -
		       ($bandwidth[$_]/$total_bandwidth))**2 for (0..$num_routers-1);
		print "$t\n";
	}
}

if($print_metrics)
{
	print "@{$observed[$_]}\n" for (0 .. $#observed);
}

sub fairness
{
	my ($rroutes, $rflows, $rbandwidth) = @_;
	my @routes = @$rroutes; # the routers each flow passes through
	my @flows = @$rflows;	# contains a list of flows passing through each
							# router
	my @bandwidth = @$rbandwidth; 	# contains the remaining bandwidth at
									# each router
	my @flow_bws; #how much bandwidth each flow gets
	#print "@{$flows->[$_]}\n" for (0 .. @bandwidth-1);
	my @order = ( 0 .. $#bandwidth ); # original router 
	my $total = 0;
	while(@bandwidth)
	{
		#get rid of routers with no active flows passing through them 
		for(my $router = $#bandwidth; $router >= 0; $router--)
		{
			unless(defined $flows[$router] and @{$flows[$router]})
			{
				#No flows for router $order[$router], deleting 
				splice @bandwidth, $router, 1;
				splice @flows, $router, 1;
				splice @order, $router, 1;
			}
		}
		next unless @bandwidth;

		my $min = 0;
		for my $i (1 .. $#bandwidth)
		{
			if(($bandwidth[$i]/@{$flows[$i]}) <
				($bandwidth[$min]/@{$flows[$min]}))
			{
				$min = $i;
			}
		}

		# $min now holds the index of the bottleneck router; increment all
		# remaining, active flows by the bottleneck amount ...
		for my $flow (0 .. ($num_flows*$num_paths) - 1)
		{
			if(defined $routes[$flow] )
			{
				$flow_bws[$flow] += $bandwidth[$min]/@{$flows[$min]};
			}
		}

		# ... and decrement all routers by the bottleneck amount for each active
		# flow passing through them
		my $increment = $bandwidth[$min]/@{$flows[$min]};
		for my $i (0 .. $#bandwidth)
		{
			
			$bandwidth[$i] -= $increment*@{$flows[$i]};
			$bandwidth_used[$order[$i]] += $increment*@{$flows[$i]};
			if($bandwidth[$i] < -1)
			{
				print "Router $i overcommitted ($bandwidth[$i])!\n";
			}
		}
		
		# which flows go via the bottleneck router?
		for my $flow (@{$flows[$min]})
		{
			#Flow $flow now maxed out; remove it ...
			foreach my $host (@{$routes[$flow]})
			{
				my $host_i= 0;
				$host_i++ while $order[$host_i] != $host;
				#... from $host 
				if($host_i!= $min)
				{
					# remove from the flow lists of the rest of the path
					@{$flows[$host_i]} = grep {$_ != $flow} @{$flows[$host_i]};
				}
			}
			$routes[$flow] = undef;
		}
		$flows[$min] = [];
	}
	for my $flow (0 .. $#flow_bws)
	{
		my @path = @{@$rroutes[$flow]};
		#print "@path\n";
		for my $i (0 .. $route_len - 2)
		{
			$just_observed[$path[$i]][$path[$i]] +=$flow_bws[$flow];
			$just_observed[$path[$i]][$path[$i+1]] +=$flow_bws[$flow];
			$just_observed[$path[$i+1]][$path[$i]] +=$flow_bws[$flow];
		}
		$just_observed[$path[$route_len-1]][$path[$route_len-1]]
													+=$flow_bws[$flow];

	}
	$total+=$bandwidth_used[$_] for (0 .. $num_routers-1);
	return ($total, \@flow_bws);
}

sub throughput
{
	my ($flows, $objective, @bandwidth) = @_;
	# Start the solver
	#my $lp = Math::GLPK->new(scalar @bandwidth, scalar @$objective);
	#$lp->prob_name("aggregate bandwidth maximizer");
	#$lp->obj_dir($LPX_MAX);

	# set all variables to be non-negative
	#$lp->set_all_var_bnds($LPX_LO, 0);

	# defined and add the objective (row vector)
	#$lp->objective( Math::MatrixSparse->new_from_cols( $objective ) );

	# define and add the right-hand side (column vector)
	#$lp->rhs( Math::MatrixSparse->new_from_rows( \@bandwidth ) );

	# define and add the constraint matrix
	# do this from cols so that each col is one flow -- that way we can pick
	# them one at a time and just push them onto the list
	#$lp->constraint( Math::MatrixSparse->new_from_cols( $flows ) );

	# use the presolver
	#$lp->int_parm($LPX_K_PRESOL, 1);

	# apply scaling, before solving the problem
	#$lp->int_parm($LPX_K_SCALE, 1);
	#$lp->int_parm($LPX_K_MSGLEV, 0);

	# solve the problem by the simplex method
	#$lp->simplex;

	# check if optimal solution is obtained
	#$lp->is_optimal or die "Could not obtain optimal solution!\n";

	#my @results;
	#my $i = 0;
	#my $total = 0;
	#foreach (split "\n", $lp->get_opt_sol->as_pretty_string)
	#{
	#	s/[^.\d]//g;
	#	$results[$i++] = $_ + 0;
	#	last if $i >= ($num_flows * $num_paths);
	#}
	#foreach my $flow (0 .. ($num_flows * $num_paths) -1)
	#{
	#	my $last = undef;
	#	foreach my $router (0 .. $#bandwidth)
	#	{
	#		if($flows->[$flow]->[$router])
	#		{
	#			$bandwidth_used[$router] += $results[$flow] ;
	#			$total += $results[$flow];
	#			$just_observed[$router][$router] +=$results[$flow];
	#			if(defined $last)
	#			{
	#				$just_observed[$router][$last] += $results[$flow];
	#				$just_observed[$last][$router] += $results[$flow];
	#			}
	#			$last = $router;
	#		}
	#	}
	#}#
	#return ($total, \@results);
}

sub get_route_col
{
	my @col = ();
	push @col, 0 for (0 .. $#bandwidth);
	$col[$_] = 1 for (&get_route);
	return @col;
}
sub set_selection_level
{
	#All paths associated with the same flow get the same selection_level
	return if ++$cur_path < $num_paths;
	$cur_path = 0; # OK, new flow, reset path counter
	if($weight_rand_level)
	{
		my $r = rand;
		my $t = $weights[0];
		$s = 0;
		while ($t <= $r)
		{
			$s++;
			$t+=$weights[$s];
		}
	}
	elsif($rand_level)
	{
		$s = int(rand(16));
	}
}

sub get_route
{
	my %route;
	$guard_selection = 1;
	&update_metrics() if $proposed_metric;
	&set_selection_level() if $weight_rand_level or $rand_level;
	while(keys %route < $route_len)
	{
		$route{&get_router()} = 1;
		$guard_selection = 0;
	}
	return keys %route;
}

sub get_router
{
	# no data yet, just choose randomly
	unless($total_metric > 0)
	{
		return int(rand($num_routers));
	}

	my $router = 0;

	#flat random (no qualifications for entry positions)
	if($flat_selection)
	{
		return int(rand($num_routers));
	}

	#current Tor
	if($tor_selection)
	{
		if($unrated > 0.9*$num_routers)
		{
			#this is really only the case if we're using the proposed
			#metric, beause it take a long time to bootstrap...
			return int(rand($num_routers));
		}
		do
		{
			$router = 0;
			my $bw = rand($total_metric);
			$bw -= $metrics[$router++] while($bw > $metrics[$router]);
		} while ($guard_selection and ($metrics[$router] < $guard_thresh));
		return $router;
	}

	#proposed tor
	if($proposed_selection)
	{
		if($s == 0)
		{
			my $min = $guard_selection ? $num_routers/2 : 0;
			$router = $min + int(rand($num_routers - $min));
		}
		elsif(not $guard_selection and (rand($num_routers) < $unrated))
		{
			#pick randomly from the unrated routers to bootstrap them
			$router = $#sorted_indices-int(rand($unrated));
		}
		else
		{
			#otherwise, pick from the rated routers by selection level
			do
			{ $router=int(($num_routers-$unrated)*(2**($s*rand())-1)/(2**$s-1));
			} while ($guard_selection and ($router > $num_routers/2));
		}
		return $sorted_indices[$router];
	}
}

sub update_metrics
{

	#Unachievable, ideal case
	if($bandwidth_metric)
	{
		@metrics = @bandwidth;
	}

	#current case
	if($tor_metric)
	{
		@metrics = map {$BW_MAX < $observed[$_][$_]  ?
						$BW_MAX :
						$observed[$_][$_]
						} (0..$num_routers - 1);
	}

	if($proposed_metric)
	{
		my %rand_peers;
		while(keys %rand_peers < $N_PEERS)
		{
			$rand_peers{int(rand($num_routers))}=1;
		}
		@metrics = map {&median(@{$_}[keys %rand_peers])} @observed;
	}

	if($eigen_metric)
	{
		my @oldmetrics;
		my $delta;
		my @row_totals = 0 for (0..$num_routers-1);
		for my $i (0..$num_routers-1)
		{
			$row_totals[$i] += $observed[$i][$_] for (0..$num_routers-1);
			$row_totals[$i] = 1 unless $row_totals[$i];
		}

		unless($total_metric)
		{
			$oldmetrics[$_] = 1/$num_routers for (0..$num_routers-1);
			print STDERR "Starting EigenSpeed from uniform vector...\n" if $debug;
		}
		else
		{
			@oldmetrics=@metrics;
		}
		my $iter=0;
		do
		{
			for my $i (0..$#oldmetrics)
			{
				$metrics[$i] = 0;
				$metrics[$i] += $observed[$_][$i]*$oldmetrics[$_]/
								$row_totals[$_] for (0..$#oldmetrics);
			}
			$delta = 0;
			$delta += ($metrics[$_]-$oldmetrics[$_])**2 for (0..$#metrics);
			@oldmetrics = @metrics;
			$iter++;
		} while($delta > 1e-10 and $iter != $max_eigen_iter);
		print STDERR "EigenSpeed took $iter iterations.\n" if $debug;
	}

	@sorted_indices = sort {$metrics[$b]<=>$metrics[$a]} (0..$num_routers-1);
	$guard_thresh = $metrics[$sorted_indices[$num_routers/2]];
	$total_metric = 0;
	$total_metric += $_ for @metrics;
	$unrated = 0;
	$unrated++ while @sorted_indices > $unrated and
			$metrics[$sorted_indices[$#sorted_indices - $unrated]] == 0;
	print STDERR "Unrated: $unrated\n" if $debug;
}

sub median
{
	my $t=0;
	for (@_)
	{
		$t++ unless $_;
	}
	return 0 if $t > 2;
	for(my $i=0; $i < $N_PEERS; $i++)
	{
		$t=0;
		for(my $j=0; $j < $N_PEERS; $j++)
		{
			next if $i == $j;
			$t += $_[$i] <=> $_[$j];
		}
		return $_[$i] unless $t;
	}
	#something is wrong, fall back on known good
	print STDERR "first-pass median failed, using sort (@_)\n" if $debug;
	return (sort {$a <=> $b} @_)[($N_PEERS-1)/2];
}

sub combine
{
	my ($old, $new) = @_; 
	return ($new > $update_thresh * $old ? $new : $old) if $aggr_max;
	if($aggr_ewma)
	{
		return($old ? ($alpha*$old +(1-$alpha)*$new) : $new);
	}
	if($aggr_mwma)
	{
		#Like an EWMA, only the larger value gets the higher weight...
		return($old ? ($old+$new)/2 + (0.5-$alpha)*abs($old-$new) : $new);
	}
}

sub node_crash($)
{
	my $uptime = $_[0];
	return (rand(1) < $crash_prob);
}
sub min(@)
{
	my @vals = @_;
	my $min = 2**32;
	$min = $min < $_ ? $min : $_ for (@vals);
	return $min;
}
