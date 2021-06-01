open(FA, $ARGV[0]);
while (<FA>) {
	chomp();
	@arr = split("\t", $_);
	if($.==1){
		next;
	}
	# print $arr[1],"\|\t";
	# $arr[1] =~ s/^\s+//;
	# $arr[1] =~ s/\s+$//;
	# print $arr[1],"\|\t";
	$sampleName = "";
	if($arr[1] =~ /(.*)\:(\d+)$/){
		$sampleName = $1;
		# print $sampleName,"\|\t";
		$sampleName =~ s/^\s+//;
		$sampleName =~ s/\s+$//;
		# print $sampleName,"\|\n";
		$samples{$sampleName}++;
	}
}
close FA;

open(FA, $ARGV[0]);
while (<FA>) {
	chomp();
	@arr = split("\t", $_);
	if($.==1){
		print $_,"\tAvailability\n";
		next;
	}
	if($arr[0] eq ""){
		next;
	}
	for($i=0;$i<=$#arr;$i++){
		$arr[$i] =~ s/^\s+//;
		$arr[$i] =~ s/\s+$//;
	}
	if(exists $samples{$arr[3]} and exists $samples{$arr[4]}){
		print $_,"\t","P1P2","\n";
		next;
	}
	if(exists $samples{$arr[3]}){
		print $_,"\t","P1","\n";
		next;
	}
	if(exists $samples{$arr[3]}){
		print $_,"\t","P2","\n";
		next;
	}
	print $_,"\t","NA","\n";
}
close FA;