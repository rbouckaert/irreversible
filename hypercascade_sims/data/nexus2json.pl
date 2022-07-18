open(FIN,"truth.nwk") or die "Cannot open file truth.nwk for writing";
<FIN>;
for ($i = 0; $i < 100; $i++) {
	$s = <FIN>;
	chomp($s);
	$tree[$i] = $s;
}
close FIN;

for ($i = 0; $i < 100; $i++) {
	open(FOUT,">data$i.json") or die "Cannot open file data$i.json for writing";
	print FOUT "{\n";
	print FOUT "run:\"$i\",\n";
	print FOUT "tree:\"$tree[$i]\",\n";
	print FOUT "sequences:\"\n";
	process($i);
	print FOUT "\"\n";
	print FOUT "}\n";
	close FOUT;
}


sub process {
	open(FIN, "Tree".($i+2)."_Barcodes.nex");
	while ($s = <FIN>) {
		if ($s =~ /(taxon[0-9]+)\s+(.*)\r/) {
			$taxon = $1;
			$sequence = $2;
			print FOUT "<sequence taxon='$taxon' totalcount='4' value='$sequence'/>\n";
		}
	}
	close FIN;
}