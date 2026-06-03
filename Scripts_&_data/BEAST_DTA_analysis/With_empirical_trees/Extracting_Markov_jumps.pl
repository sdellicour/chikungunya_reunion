#!/usr/bin/perl

$num_args = $#ARGV + 1;
if (($num_args != 1) && ($num_args != 3)) {
	print "Usage: collect_times <burning #> [from_state] [to_state]\n";	
	exit(-1);
}

$burn_in = $ARGV[0];
$print_all = 1;
if ($num_args == 3) {
	$print_all = 0;
	$from = $ARGV[1];
	$to = $ARGV[2];
} else {
	$from = "ALL";
	$to = "ALL";
}

print STDERR "Removing states < $burn_in as burn in.\n";
print STDERR "Recording jumps from $from to $to.\n";

for ($i = 0; $i < 3; $i++) {
	$line = <STDIN>;
}

my $stateCount=0;
my $jumpCount=0;

print "state\tfrom\tto\ttime\n";
while($line = <STDIN>) {
	chomp($line);
	($state, $count, $jumpf) = split('\t',$line,3);
	if ($state >= $burn_in) {
		$stateCount++;
		$jumpf =~ s/\s\d+$//; #There is a redundant total count at the end of the jump string that we need to remove
#		$jumpf =~ s/\{//;
#		$jumpf =~ s/\}//;
#		print STDERR "$jumpf\n";
		@jumps = split('},{',$jumpf);
		$n_jumps = scalar(@jumps);
		for ($i = 0; $i < $n_jumps; $i++) {
			$tmp = $jumps[$i];
			$tmp =~ s/\{\{//;
			$tmp =~ s/\}\}//;
#			print STDERR "$tmp\n";
			($site, $time, $ori, $dest) = split(',',$tmp,4);
			#$info[2] =~ s/\s\d+$//;
			if (
				($print_all == 1) ||
				(($ori eq $from) && ($dest eq $to))
			) {
				$jumpCount++;
				print "$state\t$ori\t$dest\t$time\n";
			}
		}
	}	
}
print STDERR "total state counts = $stateCount , total jumps = $jumpCount\n";
