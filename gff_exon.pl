#!/usr/bin/perl
# (1) Define and Initialize variables.
# (2) Help menu.
# (3) Hash of length of all transcripts. 
# (4) Hash of total exon length.
# (5) Hash of all GFF features. 
# (6) Switch for different functions.
# (7) Functions.
#	(7.1) Hash of length of all transcripts.
#		Open Gff File and create a Hash of Arrays. The first key of the hash
#		corresponds to the ID of the gene the second key to the transcript ID
#		and the values are an array of the transcript lengths.
#	(7.2) Hash of total exon length.
#		Transverse the hash and create a Hash of Arrays. The first key of the hash
#		corresponds to the ID of the gene the second key to the transcript ID
#		and the values are the sum total of elements in the array (total exon
#		length).
#	(7.3) Hash of all GFF features.
#		Open GFF File and create a Hash of Arrays. For this hash the first key
#		corresponds to the ID of the transcripts and the values of the array are
#		the GFF lines corresponding to said tranacript.
#	(7.4) Shortest transcript features.
#		Traverse the Hash of total exon length and determine the shortest exon  
#		transcript for each gene. Once found, traverse the Hash of all features  
#		and retrieve all features corresponding to said transcript.
#	(7.5) Longest transcript features.
#		Traverse the Hash of total exon length and determine the longest exon  
#		transcript for each gene. Once found, traverse the Hash of all features  
#		and retrieve all features corresponding to said transcript.
#	(7.6) Sort transcripts in descending order.
#		Traverse the Hash of total exon length and sort the secondary keys in 
#		descending order.
#	(7.7) Sort transcripts in ascending order.
#		Traverse the Hash of total exon length and sort the secondary keys in 
#		ascending order.  
#	(7.8) Shortest Gene transcript.
# 		Traverse the Hash of total exon length and find the shortest exon in 	
#		the whole GFF file. Output gene and transcript ID. Option to print 
#		all GFF features associated to transcript depending command line 
#		response. Traverse Hash of all GFF features to retrieve information.
#	(7.9) Longest Gene transcript.
# 		Traverse the Hash of total exon length and find the longest exon in 	
#		the whole GFF file. Output gene and transcript ID. Option to print 
#		all GFF features associated to transcript depending command line 
#		response. Traverse Hash of all GFF features to retrieve information.
#	(7.10) Statistical analysis.
#		Carry out statistical analysis of all exons in GFF file to calculate
#		number of total exons, mean, standard deviation and variance. Traverse
#		Hash of total exons and push secondary key values into an array. 
#	(7.11) Plot differences of longest and shortest exon for each gene.
#		Traverse hash of exons and calculate shortest and longest exon for
#		each. Compute difference and represent value in logarithmic scale
#		base 2. Only output genes with different exon length. 
#	(7.12) Search for gene.
#		Search for gene of interest by traversing Hash of exons and outputting
#		exon length. Option to print all GFF features associated to transcript 
#		depending command line response. Traverse Hash of all GFF features to 
#		retrieve information. This works because transcript ID starts with 
#		gene ID. 

# (1)
use strict;
use warnings;
my ($OPTION,$infile, $gene_id) = @ARGV;
my ($gff_function);
# (2)
if (!defined $infile or $OPTION =~ m/-h/ or $OPTION !~ m/-gff_s/ && $OPTION !~ m/-gff_l/ && $OPTION !~ m/-sort_a/ && $OPTION !~ m/-sort_d/ && $OPTION !~ m/-exon_s/ && $OPTION !~ m/-exon_l/ && $OPTION !~ m/-stat/ && $OPTION !~ m/-plot/ && $OPTION !~ m/-search/){ print "Usage: ./gff_length.pl [OPTION] [filename.gff]\n";
	            print "-h\n\t Help\n";
	            print "-gff_s\n\t Output GFF of SHORTEST transcripts of each gene\n";	 
		    print "-gff_l\n\t Output GFF of LONGEST transcripts of each gene\n";
		    print "-sort_d\n\t Print genes and exon lengths sorted in DESCENDING order\n";
		    print "-sort_a\n\t Print genes and exon lengths sorted in ASCENDING order\n";
		    print "-exon_s\n\t Output Gene ID Transcript ID of the SHORTEST EXON of all the file\n";
		    print "-exon_l\n\t Output Gene ID Transcript ID of the LONGEST EXON of all the file\n";
		    print "-stat\n\t Output statistical analysis of all exon lengths\n";
		    print "-plot\n\t Output HISTOGRAM in logarithmic scale of difference between LONGEST and SHORTEST exon for each gene - genes with no difference excluded\n";
		    print "-search\n\t search GENE ID and output EXON lengths with additional OPTION of outputting all associated gff features\n\t Usage: ./gff_length.pl -search [filename.gff] [Gene_ID] \n";
exit(0);
};
# (3)
my (%Hash_transcripts)=read_file($infile);
# (4)
my(%Hash_exons)=exon_extractor(%Hash_transcripts);
# (5)
my (%Hash_all_features)=all_feature_extractor($infile);
#(6)
SWITCH: {
    ($OPTION =~ /-gff_s/io) && do {
        $gff_function = \&sort_shortest_gff;
        last SWITCH;
    };
    ($OPTION =~ /-gff_l/io) && do {
        $gff_function = \&sort_longest_gff;
        last SWITCH;
    };
    ($OPTION =~ /-sort_a/io) && do {
        $gff_function = \&sort_ascending;
        last SWITCH;
    };
    ($OPTION =~ /-sort_d/io) && do {
        $gff_function = \&sort_descending;
        last SWITCH;
    };
    ($OPTION =~ /-exon_s/io) && do {
        $gff_function = \&shortest_exon;
        last SWITCH;
    };
    ($OPTION =~ /-exon_l/io) && do {
        $gff_function = \&longest_exon;
        last SWITCH;
    };    
    ($OPTION =~ /-stat/io) && do {
        $gff_function = \&statistics;
        last SWITCH;
    };        
    ($OPTION =~ /-plot/io) && do {
        $gff_function = \&plot;
        last SWITCH;
    };
    ($OPTION =~ /-search/io) && do {
        $gff_function = \&grep_exon_feature;
        last SWITCH;
    };            
};
&$gff_function(\%Hash_exons, \%Hash_all_features, $gene_id);;
exit 0;
# (7)
#(7.1)
sub read_file {
my %hash_arrays;
my $input = shift;
open(FILE_A, $input) || die ("Couldn't open File A");
while (<FILE_A>){
		next if /^\#/o;
		my @line = split (' ',$_);
		$line[9] =~ s/^"//g;
		$line[9] =~ s/";$//g;
		$line[11] =~ s/^"//g;
		$line[11] =~ s/"$//g;
		if($line[2] =~ m/utr|First|Internal|Terminal|Single/g){push @{$hash_arrays{$line[9]}{$line[11]}}, $line[4]-$line[3] + 1};
	
	};
	
close (FILE_A);
return %hash_arrays;
}
#(7.2)
sub exon_extractor {
my (%old_hash)= @_;
my %new_hash;
foreach my $key ( keys %old_hash ) {
	foreach my $key2 (keys %{ $old_hash{$key} }){
		my $total_length = 0;
		my @features = ();
		my $string1 = ();            
			foreach ( @{$old_hash{$key}{$key2}})  {                   
				$total_length += $_ ;
	
			}
			foreach ( @{$old_hash{$key}{$key2}}[1] )  {                   
				push @features, $_;               
			}
push@{$new_hash{$key}{$key2}}, $total_length;
	}
}
return %new_hash;
}
#(7.3)
sub all_feature_extractor{
my $input = shift;
my %hash_arrays;
open(FILE_A, $input) || die ("Couldn't open File A");
while (<FILE_A>){
		next if /^\#/o;
		my @line = split (' ',$_);
		$line[11] =~ s/^"//g;
		$line[11] =~ s/"$//g;
		if($line[2] =~ m/utr|First|Internal|Terminal|Single/g){
			push @{$hash_arrays{$line[11]}}, $_;
		}
	
	};
	
	
	
close (FILE_A);
return %hash_arrays;
}
# (7.4)
sub sort_shortest_gff(\%%) {
my %Hash_1 = %{ $_[0] }; my %Hash_2 = %{ $_[1] };  my $grep_id = $_[2];
foreach my $oh (keys %Hash_1) {
    my $lowest = 99999999999;
    my $lowest_key;
    my $value;
    foreach my $ih (keys %{$Hash_1{$oh}}) {
        $value = @{$Hash_1{$oh}{$ih}}[0];
        if ($value < $lowest) {
           	$lowest = $value;
            	$lowest_key = $ih;
	}
	}	
print "# Gene: $oh\tTranscript: $lowest_key\tExon_Length: $lowest\n";

foreach my $oh (keys %Hash_2) {
		if ($oh =~ m/$lowest_key/) {
		print $lowest."\t";
		print @{$Hash_2{$oh}};	
		}
	}
}
}
# (7.5)
sub sort_longest_gff(\%%) {
my %Hash_1 = %{ $_[0] }; my %Hash_2 = %{ $_[1] };  my $grep_id = $_[2];
foreach my $oh (keys %Hash_1) {
my $large_val = 0;
my $large_key;
my $value;
foreach my $ih (keys %{$Hash_1{$oh}}) {
	$value = @{$Hash_1{$oh}{$ih}}[0];
        if ($value > $large_val) {
           	$large_val = $value;
            	$large_key = $ih;
	}
	}	
print "# Gene: $oh\tTranscript: $large_key\tExon_Length: $large_val\n";

foreach my $oh (keys %Hash_2) {
		if ($oh =~ m/$large_key/) {
		print $large_val."\t";
		print @{$Hash_2{$oh}};
		}
	}
}
}
# (7.6)
sub sort_descending {
my %hash = %{ $_[0] }; my %hash2 = %{ $_[1] };  my $grep_id = $_[2];
foreach my $oh (keys %hash) {
    foreach my $ih (sort { @{$hash{$oh}{$b}}[0] <=> @{$hash{$oh}{$a}}[0] } keys %{$hash{$oh}}) {
        print join(", ", $oh, $ih, @{$hash{$oh}{$ih}}[0]) . "\n";
    }
}
}

# (7.7)
sub sort_ascending {
my %hash = %{ $_[0] }; my %hash2 = %{ $_[1] };  my $grep_id = $_[2];
foreach my $oh (keys %hash) {
    foreach my $ih (sort { @{$hash{$oh}{$a}}[0] <=> @{$hash{$oh}{$b}}[0] } keys %{$hash{$oh}}) {
        print join(", ", $oh, $ih, @{$hash{$oh}{$ih}}[0]) . "\n";
    }
}
}
# (7.8)
sub shortest_exon(\%%) {
my %Hash_1 = %{ $_[0] }; my %Hash_2 = %{ $_[1] }; my $grep_id = $_[2];
my $large_val = 99999999999;
my ($large_key, $value, $gene);
foreach my $oh (keys %Hash_1) {   
    foreach my $ih (keys %{$Hash_1{$oh}}) {
        $value = @{$Hash_1{$oh}{$ih}}[0];
        if ($value < $large_val) {
        	$gene = $oh;
           	$large_val = $value;
            	$large_key = $ih;
	}
	}	
}
print "# Gene: $gene\tTranscript: $large_key\tExon_Length: $large_val\n";
print "Do you wish to see all the associated gff features to this gene? YES or NO.\n";
    while (<STDIN>) {
	chomp;
	if($_ eq "YES"){
		foreach my $oh (keys %Hash_2) {
		if ($oh =~ m/$large_key/) {
		print @{$Hash_2{$oh}};
		}
	}
		exit 0;	
	}elsif($_ eq "NO"){
		exit 0;
	}else{print STDOUT "Please type YES or NO.\n";};
	
}

}
# (7.9)
sub longest_exon(\%%) {
my %Hash_1 = %{ $_[0] }; my %Hash_2 = %{ $_[1] };  my $grep_id = $_[2];
my $large_val=0;
my $large_key;
my $value;
my $gene;
foreach my $oh (keys %Hash_1) {   
    foreach my $ih (keys %{$Hash_1{$oh}}) {
        $value = @{$Hash_1{$oh}{$ih}}[0];
        if ($value > $large_val) {
           	$gene = $oh;
           	$large_val = $value;
            	$large_key = $ih;
	}
	}	
}
print "# Gene: $gene\tTranscript: $large_key\tExon_Length: $large_val\n";
print "Do you wish to see all the associated gff features to this gene? YES or NO.\n";
    while (<STDIN>) {
	chomp;
	if($_ eq "YES"){
		foreach my $oh (keys %Hash_2) {
		if ($oh =~ m/$large_key/) {
		print @{$Hash_2{$oh}};
		}
	}
		exit 0;	
	}elsif($_ eq "NO"){
		exit 0;
	}else{print STDOUT "Please type YES or NO.\n";};
	
}
}
# (7.10)
sub statistics(\%%) {
my %Hash_1 = %{ $_[0] }; my %Hash_2 = %{ $_[1] }; my $grep_id = $_[2];
my (@array, $sum, $total_diff);
foreach my $oh (keys %Hash_1) {   
    foreach my $ih (keys %{$Hash_1{$oh}}) {
	push @array, @{$Hash_1{$oh}{$ih}}[0];
	}
}	
foreach (@array) {
        $sum += $_;
    }
my $mean = ($sum/ @array);
my $n = @array;
#print $mean $n;
foreach my $exon (@array){
  my $diff = $exon - $mean;
  $diff *= $diff;
  $total_diff += $diff;
}
my $var = ($total_diff / ($n-1));
my $sd = sqrt($var);
print " ################## Statistical analysis of all exons ##################\n";
printf("Exon_count: %.d\t Mean: %.2f\t Standar_deviation: %.2f\t Variance: %.2f\n", $n, $mean, $sd, $var);
}
# (7.11)
sub plot {
my %Hash_1 = %{ $_[0] }; my %Hash_2 = %{ $_[1] }; my $grep_id = $_[2];
foreach my $oh (keys %Hash_1) {
my $large_val = 0;
my ($large_key, $value);
foreach my $ih (keys %{$Hash_1{$oh}}) {
	$value = @{$Hash_1{$oh}{$ih}}[0];
        	if ($value > $large_val) {
           		$large_val = $value;
            		$large_key = $ih;
		}
	}
	my $lowest = 99999999999;
	my $lowest_key;	
	foreach my $ih (keys %{$Hash_1{$oh}}) {
       		$value = @{$Hash_1{$oh}{$ih}}[0];
       	 	if ($value < $lowest) {
           		$lowest = $value;
            		$lowest_key = $ih;
		}
	}
my $n = ($large_val - $lowest);
if ($n != 0) {print $oh."\t"; print '▦' x (log($n)/ log(2))."\n";};
}
}
# (7.12)
sub grep_exon_feature (\%$) {
my %Hash_1 = %{ $_[0] }; 
my %Hash_2 = %{ $_[1] }; 
my $grep_id = $_[2];
foreach my $oh (keys %Hash_1) {
if ($oh =~ m/$grep_id/) {
	print "Gene:$grep_id\tExon lengths: ";
    foreach my $ih (sort { @{$Hash_1{$oh}{$b}}[0] <=> @{$Hash_1{$oh}{$a}}[0] } keys %{$Hash_1{$oh}}) {
 	print @{$Hash_1{$oh}{$ih}}[0]."  ";   
    }
    print "\nDo you wish to see all the associated gff features to this gene? YES or NO.\n";
    while (<STDIN>) {
	chomp;
	if($_ eq "YES"){
		foreach my $oh (keys %Hash_2) {
		if ($oh =~ m/$grep_id/) {
		print @{$Hash_2{$oh}};
		}
	}
		exit 0;	
	}elsif($_ eq "NO"){
		exit 0;
	}else{print STDOUT "Please type YES or NO.\n";};
	
}
}
}
	print "ERROR $grep_id not found!\n";
	exit 0; 	
}
