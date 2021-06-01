open(FA, $ARGV[0]);
while(<FA>){
	chomp();
	@arr = split("\t", $_);
	if($arr[2] ne "" and $arr[3] ne ""){
		$parents{$arr[2]}++;
		$parents{$arr[3]}++;	
	}
	if($arr[1]=~/(.*)\:(\d+)/){
		$gt = $1;
		$genotypes{$gt}++;
	}
	
}
close FA;


open(FA, $ARGV[0]);
while(<FA>){
	chomp();
	@arr = split("\t", $_);
	if($arr[1]=~/(.*)\:(\d+)/){
		$gt = $1;
		if($parents{$gt}){
			if((exists($genotypes{$arr[2]})) and (exists($genotypes{$arr[3]}))){
				print "SPECIAL\t",$arr[0],"\t",$arr[1],"\t",$gt,"\t--\t--","\n";
			}
			else{
				print "P\t",$arr[0],"\t",$arr[1],"\t",$gt,"\t--\t--","\n";	
			}
			
		}
		else{
			if((exists($genotypes{$arr[2]})) and (exists($genotypes{$arr[3]}))){
				print "F1\t",$arr[0],"\t",$arr[1],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t","\n";	
			}
			else{
				if($arr[2] eq "" and $arr[3] eq ""){
					print "X_P\t",$arr[0],"\t",$arr[1],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t","\n";
				}
				elsif((!exists($genotypes{$arr[2]})) and (!exists($genotypes{$arr[3]}))){
					print "X_P1_P2\t",$arr[0],"\t",$arr[1],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t","\n";
				}
				elsif((!exists($genotypes{$arr[2]}))){
					print "X_P1\t",$arr[0],"\t",$arr[1],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t","\n";
				}
				elsif((!exists($genotypes{$arr[3]}))){
					print "X_P2\t",$arr[0],"\t",$arr[1],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t","\n";
				}
				else{
					print "MISLEAD\t",$arr[0],"\t",$arr[1],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\t","\n";
				}
			}
			
		}
	}
	else{
		print $_,"\n";
	}
}
close FA;