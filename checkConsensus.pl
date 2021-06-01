$IUPAC{'A'} = A;
$IUPAC{'T'} = T;
$IUPAC{'G'} = G;
$IUPAC{'C'} = C;


open(FA, $ARGV[0]);
while (<FA>) {
	chomp();
	@arr = split("\t", $_);
	if($.==1){
		@head = @arr;
		next;
	}
	for ($i=11;$i<=$#head;$i++){
		$hash{$arr[0]}{$head[$i]} = $arr[$i]
	}
}
close FA;

open(FA, $ARGV[1]);
while (<FA>) {
	chomp();
	@arr = split("\t", $_);
	$length = $#arr;
	foreach $mIDs(sort keys %hash){
		$base_a = $hash{$mIDs}{$arr[1]};
		$base_b = $hash{$mIDs}{$arr[2]};
		# if($base_a ne $base_b && exists $IUPAC{$base_a} && exists $IUPAC{$base_b}){
		if($base_a ne $base_b && $base_a ne "N" && $base_b ne "N"){			
			print $arr[0],"\t",$mIDs,"\t",$base_a,"\t",$base_b,"\n";	
		}
		
	}
	
}
close FA;