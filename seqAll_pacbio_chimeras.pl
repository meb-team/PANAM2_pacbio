use strict;
use warnings;

# USAGE
# perl seqAll.pl /dir/containing/<18S fasta file>/ amplicons_709.txt



my $in = $ARGV[0];
#~ my $out = $ARGV[1];

my $sample_filename = $ARGV[1];
my $sample_name;

my $i = 0;	# seq number
my $ab = 0;	# abundance of preclusters

# print "$in\n";

print "$sample_filename\n";

open (Fsample, "$sample_filename") or die ("error sample file name\n");
while (my $lisample = <Fsample>){
	chomp ($lisample);
	
	my @tsample = split("\t", $lisample);
	
	# $sample_name = "$tsample[1]";
	# print "$in\n\n";
	# print "$tsample[1]\n";
	if($tsample[1] eq $in){
		$sample_name = $tsample[0];
		# print "$sample_name\n";
	}

}
close Fsample;


open (SEQALL, '>', $in."/".$sample_name.".final.fasta") or die "Couldn't open: SEQALL"; 
open (CORRES, '>', $in."/".$sample_name.".final.list") or die "Couldn't open: CORRES";



print "$in"."/"."$in".".count_table\n";

open (F0, $in."/all18S.filtered.fasta") or die ("error in.fasta\n");
while (my $li0 = <F0>){
	chomp ($li0);
	
	my $premc = substr($li0,0,1);
	
	if($premc eq '>'){
		$i ++;
		
		open(F1, $in."/".$in.".count_table") or die ("error in.count_table\n");
		while (my $li1 = <F1>){
			chomp ($li1);
			
			my @tcount = split("\t", $li1);
			
			my $id = substr($li0,1,100);
			if($li1 =~ m/$id/){
					$ab = $tcount[1];
			}
		}
		close F1;
		
		print SEQALL ">".$sample_name."_Seq".$i."-".$ab."\n";
		print CORRES substr($li0,1,100)."\t"."Seq".$i."-".$ab."\n";
	}
	else { print SEQALL "$li0\n"}
}
close F0;

close SEQALL;
close CORRES;