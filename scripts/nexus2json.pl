@times = (1,30,60,120,290,500);
@tcounts = (5,10,20,40,80,160,320);





foreach $taxa (@tcounts) {

open(FIN,"truth-$taxa.trees.nwk") or die "Cannot open file truth-$taxa.trees.nwk for reading";
<FIN>;
for ($i = 0; $i < 100; $i++) {
	$s = <FIN>;
	chomp($s);
	$tree[$i] = $s;
}
close FIN;

$taxnames = '';
for ($i = 1; $i <= $taxa; $i++) {
	$taxnames .= ',taxon'.$i;
}
$taxnames = substr($taxnames,1);

	foreach $time (@times) {

		print "taxa=$taxa time=$time\n";
		
		$dir = "data_".$taxa."taxa_".$time."time";
		`mkdir $dir`;

		for ($i = 0; $i < 100; $i++) {
			open(FOUT,">$dir/data$i.json") or die "Cannot open file $dir/data$i.json for writing";
			print FOUT "{\n";
			print FOUT "run:\"$i\",\n";
			print FOUT "taxa:\"$taxa\",\n";
			print FOUT "time:\"$time\",\n";
			print FOUT "tree:\"$tree[$i]\",\n";
			print FOUT "taxnames:\"$taxnames\",\n";
			print FOUT "sequences:\"\n";
			process($i);
			print FOUT "\"\n";
			print FOUT "}\n";
			close FOUT;
		}

	}
}

sub process {
	$file = "Tree".($i+2)."_taxa".$taxa."_time".$time."_Barcodes.nex";
	open(FIN, $file) or die "Cannot open Tree file $file for reading";
	for ($j = 1; $j <= $taxa; $j++) {
		$taxon[$j] = '';
	}
	while ($s = <FIN>) {
		if ($s =~ /taxon([0-9]+)\s+(.*)\r/) {
			$t = $1;
			$sequence = $2;
			$taxon[$t] = $taxon[$t].$sequence;
		}
	}
	for ($j = 1; $j <= $taxa; $j++) {
		print FOUT "<sequence taxon='taxon$j' totalcount='4' value='$taxon[$j]'/>\n";
	}
	close FIN;
}