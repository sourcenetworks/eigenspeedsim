#!/usr/bin/perl -w

use strict;
my $thresh = 1e-10;
my $liar_thresh = 1e-5;
my @T;
my $delta;
my $debug = 0;
my $min = 0;
my $max = 0;
my $avg = 0;
my $goog = 0;
my $selfrank = 1;
my $nearclique = 0;
my $N = -1;
my $n = -1;
my @files;
my @liar;
my $recalc = 1;
my $noise = 0;
my $killed = 0;
my $doliar = 0;
my $doearlyliar = 0;
my $trustednodes = 0;

$SIG{INT} = \&kill;

while($_ = shift)
{
	if($_ eq "-min") { $min = 1; }
	elsif($_ eq "-max") { $max = 1; }
	elsif($_ eq "-avg") { $avg = 1; }
	elsif($_ eq "-goog") { $goog = shift; }
	elsif($_ eq "-liar") { $doliar = 1; }
	elsif($_ eq "-earlyliar") { $doearlyliar = 1; }
	elsif($_ eq "-nearclique") { $nearclique = 1; }
	elsif($_ eq "-self") { $selfrank = 0; }
	elsif($_ eq "-trusted") { $trustednodes = 1; }
	elsif($_ eq "-limit") { $N = shift; }
	elsif($_ eq "-d") { $debug = 1; }
	elsif($_ eq "-n") { $noise = shift; }
	elsif($_ eq "-thresh") { $thresh = shift; }
	elsif($_ eq "-L") { $liar_thresh = shift; }
	else { unshift @files, $_; }
}
@ARGV = @files if @files;

while(<>)
{
	my @r = split;
	push @T, \@r;
}
if($noise)
{
	for my $i (0..$#T)
	{
		$T[$i]->[$_] *= ((1-$noise) + 2*$noise*rand) for (0..$#T);
	}
}
unless($selfrank)
{
	$T[$_]->[$_] = 0 for (0..$#T);
}

for my $i (0..$#T)
{
	for my $j (0..$i)
	{
		if($min)
		{
			$T[$i]->[$j] = $T[$j]->[$i] =
					$T[$i]->[$j] < $T[$j]->[$i] ? $T[$i]->[$j] : $T[$j]->[$i];
		}
		elsif ($max)
		{
			$T[$i]->[$j] = $T[$j]->[$i] =
					$T[$i]->[$j] > $T[$j]->[$i] ? $T[$i]->[$j] : $T[$j]->[$i];
		}
		elsif($avg)
		{
			$T[$i]->[$j] = $T[$j]->[$i] = ($T[$i]->[$j] + $T[$j]->[$i])/2; 
		}
	}
}

#normalize the rows
for my $i (0..$#T)
{
	my $sum = 0;
	$sum += $_ for (@{$T[$i]});
	if($sum)
	{
		$T[$i]->[$_] /= $sum for (0..$#T);
	}
	#begin early liar detection...
	next unless $doearlyliar;
	my $liar = 0;
	$liar += (1/$#T-$T[$i]->[$_])**4 for (0..$#T);
	if($liar >= 0.05)
	{
		print STDERR "Node $i is a liar ($liar])!\n";
		$T[$i]->[$_] = 0 for (0..$#T);
	}
}
if($nearclique)
{
	my $clique_thresh = 1/(@T*100);
	for my $i (0..$#T)
	{
		my $sum = 0;
		for my $j (0..$#T)
		{
			if($T[$i]->[$j] < $clique_thresh)
			{
				$T[$i]->[$j] = 0;
			}
			else
			{
				$sum += $T[$i]->[$j];
			}
		}
		if($sum)
		{
			$T[$i]->[$_] /= $sum for (0..$#T);
		}
	}
}

if($goog)
{
	for my $i (0..$#T)
	{
		for my $j (0..$i-1)
		{
			$T[$i]->[$j] = $goog*$T[$i]->[$j] + (1-$goog)/@T;
			$T[$j]->[$i] = $goog*$T[$j]->[$i] + (1-$goog)/@T;
		}
		$T[$i]->[$i] = $goog*$T[$i]->[$i] + (1-$goog)/@T;
	}
}

my %seen;
my @queue;
my @cliques;
my %clique;
for my $node (0..$#T)
{
	next if defined $seen{$node};
	push @queue, $node;
	while(defined($node = shift @queue))
	{
		$clique{$node} = 1;
		for my $peer (grep {$T[$node]->[$_] != 0} (0..$#T))
		{
			push @queue, $peer unless $clique{$peer};
			$clique{$peer} = 1;
		}
	}
	print STDERR "Found a clique of size ".scalar keys(%clique)."!\n" if $debug;
	push @cliques, [keys %clique];
	$seen{$_} = 1 for keys %clique;
	%clique = ();
}


my (@oldt, @t);
#@oldt = @{$T[0]};
$"="\n";
#print "@oldt";

while($recalc)
{
	if($trustednodes)
	{
		$oldt[$_] = 0 for (0..$#T-5);
		$oldt[$_] = 1/5 for ($#T-4..$#T);
	}
	else
	{
		$oldt[$_] = 1/@T for (0..$#T);
	}
	$recalc = 0;
	$killed = 0;
	$n = 0;
	do
	{
		for my $i (0..$#oldt)
		{
			$t[$i] = 0;
			$t[$i] += $T[$_]->[$i]*$oldt[$_] for (0..$#oldt);
		}
		$delta = 0;
		$delta += ($t[$_]-$oldt[$_])**2 for (0..$#t);
		@oldt = @t;
		print STDERR "$delta\n" if $debug;
	} while($delta > $thresh && $N - $n++ != 0 && not $killed);
	print STDERR ($n-1)." iterations completed, delta = $delta\n";
	
	next unless $doliar;
	for my $i (0..$#t)
	{
		$liar[$i] = 0;
		$liar[$i] += ($T[$i]->[$_]-$t[$_])**4 for (0..$#t);
		$liar[$i] -= ($T[$i]->[$i]-$t[$i])**4;
		if($liar[$i] >= $liar_thresh)
		{
			print STDERR "Node $i is a liar ($liar[$i])!\n";
			#$T[$i]->[$_] = 0 for (0..$#t);
			$T[$i]->[$_] = $T[$_]->[$i] = 0 for (0..$#t);
			$recalc = 1;
		}
	}
	#re-renormalize the rows after removing liar nodes
	for my $i (0..$#T)
	{
		my $sum = 0;
		$sum += $_ for (@{$T[$i]});
		if($sum)
		{
			$T[$i]->[$_] /= $sum for (0..$#T);
		}
	}
}

for (0..$#t)
{
	print "$t[$_]";
	print " $liar[$_]" if $doliar;
	print "\n";
}

sub kill
{
	$killed = 1;
}
