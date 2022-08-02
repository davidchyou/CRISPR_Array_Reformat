use JSON::PP qw(decode_json);

my $spacer_file = "";
my $out = "";

my $ind = 0;
foreach(@ARGV) {

	if (@ARGV[$ind] eq '-in') {
		$spacer_file = @ARGV[$ind + 1];
		if (! (-e $spacer_file)) {
			die "cannot open file: " . $spacer_file . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-out') {
		$out = @ARGV[$ind + 1];
	}
	
	$ind++;
}

if (length($out) == 0) {
	$out = "$spacer_file.gff";
}

my $flag = doGFF($spacer_file, $out);

sub doGFF {
	my ($path, $out) = @_;
	
	open(RD, $path);
	my @arr_spacer_file=<RD>;
	close(RD);
	
	my $spacer_type; 
	chomp $arr_spacer_file[0];
	$arr_spacer_file[0] =~ s/^\s+//;
	
	my $flag = 0;
	if($arr_spacer_file[0]=~ /ORGANISM/) {
		$spacer_type="crt";
		$flag = CRTToGFF($path, $out);
	} elsif($arr_spacer_file[0]=~ /pilercr/) {
		$spacer_type="pilercr";
		$flag = pilercrToGFF($path, $out);
	} elsif($arr_spacer_file[0]=~ /########################################/) {				
		$spacer_type="CRISPRFinder";
		$flag = CFToGFF($path, $out);
	} elsif($arr_spacer_file[0]=~ /^{/) {				
		$spacer_type="CRISPRCasFinder";
		$flag = CCFinderToGFF($path, $out);
	} else {	
		print "Error: Only CRT, pilercr, CRISPRCasFinder and CRISPRfinder format is supported at the moment.";
		exit;
	}
	
	return $flag;
}

sub CRTToGFF {
	my ($path, $out) = @_;
	
	open(RD, $path);
	my @arr_spacer_file=<RD>;
	close(RD);
		
	my $orgName;
	my $crispr_index;
	my $rangeStart;
	my $rangeEnd;
	my $strand = "+";
	my $app = "CRT";
	my $xref = ".";
	
	open(WR,">$out");
	for(my $i = 0; $i < scalar(@arr_spacer_file); $i++) {
		if ($arr_spacer_file[$i]=~/ORGANISM:  (.*)/){								
			$orgName = $1; 
			chomp $orgName; 
			$orgName =~ s/\r//g;
			($orgName) = ($orgName =~ /(\S+).*/);							
		} elsif ($arr_spacer_file[$i] =~ /CRISPR ([0-9]+).*Range: ([0-9]+) - ([0-9]+)/) {
			$crispr_index =$1; chomp $crispr_index;
			$rangeStart = $2; chomp $rangeStart;
			$rangeEnd = $3; chomp $rangeEnd;
			my $crispr_id = "CRISPR$crispr_index" . "_" . $rangeStart . "_" . $rangeEnd;
			my $crispr_length = abs($rangeEnd - $rangeStart) + 1;
			
			$i += 3;
			my $spacerIndex = 0;
			my $repeatIndex = 0;
			
			my %repeats = ();
			my @lines = ();
			
			while($arr_spacer_file[$i] =~ /([0-9]+)\t\t([ACTG]+)\t([ACTG]+)\t\[ ([0-9]+), ([0-9]+) ]/) {
				$spacerIndex++;
				$repeatIndex++;
				
				my $repeat_seq = $2; chomp $repeat_seq;
				my $spacer_seq = $3; chomp $spacer_seq;
				
				my $spacer_start = int($1) + int($4);
				my $spacer_end = int($1) + int($4) + int($5) - 1;
				my $repeat_start = int($1);
				my $repeat_end = int($1) + length($repeat_seq) - 1;
				my $spacer_id = "CRISPR$crispr_index" . "_" . "SPACER$spacerIndex" . "_" . $spacer_start . "_" . $spacer_end;
				my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeatIndex" . "_" . $repeat_start . "_" . $repeat_end;
				
				$i++;
				
				my $line1 = "$orgName\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$4\t$strand\t$xref\t";
				$line1 .= "ID=$repeat_id;Name=$repeat_id;Parent=$crispr_id;Note=$repeat_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
				
				my $line2 = "$orgName\t$app\tbinding_site\t$spacer_start\t$spacer_end\t$5\t$strand\t$xref\t";
				$line2 .= "ID=$spacer_id;Name=$spacer_id;Parent=$crispr_id;Note=$spacer_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
				
				push(@lines, $line1);
				push(@lines, $line2);
				
				if (exists $repeats{$repeat_seq}) {
					$repeats{$repeat_seq}++;
				} else {
					$repeats{$repeat_seq} = 1;
				}
			}
			
			while($arr_spacer_file[$i] =~ /([0-9]+)\t\t([ACTG]+)\t*/) {
				$repeatIndex++;
				my $repeat_seq = $2; chomp $repeat_seq;
				my $repeat_start = int($1);
				my $repeat_end = int($1) + length($repeat_seq) - 1;
				my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeatIndex" . "_" . $repeat_start . "_" . $repeat_end;
				my $repeat_len = length($repeat_seq);
				
				$i++;
				
				my $line1 = "$orgName\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$repeat_len\t$strand\t$xref\t";
				$line1 .= "ID=$repeat_id;Name=$repeat_id;Parent=$crispr_id;Note=$repeat_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
				
				push(@lines, $line1);
				
				if (exists $repeats{$repeat_seq}) {
					$repeats{$repeat_seq}++;
				} else {
					$repeats{$repeat_seq} = 1;
				}
			}
			
			my $reprep = "NA";
			my $freq = 0;
			
			foreach my $rep (keys(%repeats)) {
				if ($repeats{$rep} > $freq) {
					$reprep = $rep;
					$freq = $repeats{$rep};
				}
			}
			
			my $line0 = "$orgName\t$app\trepeat_region\t$rangeStart\t$rangeEnd\t$crispr_length\t$strand\t$xref\t";
			   $line0 .= "ID=$crispr_id;Name=$crispr_id;Note=$reprep;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
			
			print WR "$line0\n";
			foreach my $line (@lines) {
				print WR "$line\n";
			}
		}
	}
	close(WR);
	return 1;
}

sub pilercrToGFF {
	my ($path, $out) = @_;
	
	open(RD, $path);
	my @arr_spacer_file=<RD>;
	close(RD);
	
	my $orgName;
	my $crispr_index;
	my $rangeStart;
	my $rangeEnd;
	my $strand = "+";
	my $app = "PILER-CR";
	my $xref = ".";
	
	open(WR,">$out");
	for(my $i = 0; $i < scalar(@arr_spacer_file); $i++) {
		if($arr_spacer_file[$i] =~ />(.*)/) {
			$orgName = $1; chomp $orgName; 
			$orgName =~ s/\r//g;
			($orgName) = ($orgName =~ /(\S+).*/);
													
			my $previous_line=$arr_spacer_file[$i - 1]; chomp $previous_line;
			$previous_line =~ s/\r//;
			$previous_line =~ s/\s+$//;
								
								
			if($previous_line !~ /Array/) {										
				next;
			}
				
			$previous_line =~ s/Array //;
			$crispr_index = $previous_line;
								
			my $line_containing_first_pos=4;
			if($arr_spacer_file[$i + 3] =~ /===/) {
				$line_containing_first_pos = 4; 
			} elsif($arr_spacer_file[$i + 4] =~ /===/) {
				$line_containing_first_pos = 5; 
			}	
																
			$i = $i + $line_containing_first_pos;
			
			my $spacerIndex = 0;
			my $repeatIndex = 0;
			
			my @lines = ();
			my @seqs = ();
			my $reprep = "";
			
			while ($arr_spacer_file[$i] =~ / +([0-9]+) +([0-9]+) +[0-9\.]+ +([0-9]+) +[ACGT]+ +[-ACGT\.]+ + ([ACGT]+)/) {
				$repeat_start = $1; chomp $repeat_start;
				$repeat_end = $repeat_start + $2 - 1;
				
				$spacer_start = $1 + $2;
				$spacer_end = $spacer_start + $3 - 1;
				
				$spacer_seq = $4;
				$repeat_seq = "";
				
				$spacer_len = 
				
				$spacerIndex++;
				$repeatIndex++;
				
				if ($spacerIndex < 2) {
					$rangeStart = $repeat_start;
				}
				
				my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeatIndex" . "_" . $repeat_start . "_" . $repeat_end;
				my $spacer_id = "CRISPR$crispr_index" . "_" . "SPACER$spacerIndex" . "_" . $spacer_start . "_" . $spacer_end;
				
				my $line1 = "$orgName\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$2\t$strand\t$xref\t";
				$line1 .= "ID=$repeat_id;Name=$repeat_id;";
				
				my $line2 = "$orgName\t$app\tbinding_site\t$spacer_start\t$spacer_end\t$3\t$strand\t$xref\t";
				$line2 .= "ID=$spacer_id;Name=$spacer_id;";
				
				push(@seqs, $repeat_seq);
				push(@seqs, $spacer_seq);
				
				push(@lines, $line1);
				push(@lines, $line2);
				
				if($arr_spacer_file[$i + 2]=~ /==/) {
					if ($arr_spacer_file[$i + 1] =~ / +([0-9]+) +([0-9]+) +[0-9\.]+ +[ACGT]+ +[-ACGT\.]+ + ([ACGT]+)/) {
						my $repeat_start = $1; chomp $repeat_start;
						my $repeat_end = $repeat_start + $2 - 1;
						my $repeat_seq = "";
						
						$repeatIndex++;
						my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeatIndex" . "_" . $repeat_start . "_" . $repeat_end;
						
						my $line1 = "$orgName\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$2\t$strand\t$xref\t";
						$line1 .= "ID=$repeat_id;Name=$repeat_id;";
						
						$rangeEnd = $repeat_end;
						push(@seqs, $repeat_seq);
						push(@lines, $line1);
					}
					
					if ($arr_spacer_file[$i + 3] =~ / +([0-9]+) +([0-9]+) +([0-9]+) +([ACGT]+)/) {
						$reprep = $4;
					}
					
					last;
				}
				
				$i++;	
			}
			my $crispr_id = "CRISPR$crispr_index" . "_" . $rangeStart . "_" . $rangeEnd;
			my $crispr_length = abs($rangeEnd - $rangeStart) + 1; 
			
			my $line0 = "$orgName\t$app\trepeat_region\t$rangeStart\t$rangeEnd\t$crispr_length\t$strand\t$xref\t";
			   $line0 .= "ID=$crispr_id;Name=$crispr_id;Note=$reprep;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
			
			print WR "$line0\n";
			for (my $j = 0; $j < scalar(@lines); $j++) {
				my $seq = @seqs[$j];
				
				if ($seq eq "") {
					$seq = $reprep;
				}
				
				my $line = @lines[$j];
				$line .= "Parent=$crispr_id;Note=$seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
				print WR "$line\n";
			}	
		}
	}
	close(WR);
	return 1;
}

sub CFToGFF {
	my ($path, $out) = @_;
	
	open(RD, $path);
	my @arr_spacer_file=<RD>;
	close(RD);
	
	my $orgName;
	my $crispr_index;
	my $rangeStart;
	my $rangeEnd;
	my $strand = "+";
	my $app = "CRISPRFinder";
	my $xref = ".";
	
	open(WR,">$out");
	for(my $i = 0; $i < scalar(@arr_spacer_file); $i++) {
		my $line = $arr_spacer_file[$i]; chomp $line;
		
		if ($line =~ /^# Sequence:/) {
			$line=~ s/^# Sequence://;
								
			my ($orgName, $crispr_index, $rangeStart, $rangeEnd, $reprep);
								
			$orgName = $line; 
			$orgName =~ s/\r//g; 
			$orgName =~ s/>//;
			$orgName =~ s/\s+$//; 
			$orgName =~ s/^\s+//;
			
			$crispr_index = $arr_spacer_file[$i + 6]; chomp $crispr_index; 
			$crispr_index =~ s/# Crispr Rank in the sequence: //; 
			$crispr_index =~ s/\r//;
			
			my $crispr_start_stop_line = $arr_spacer_file[$i + 7]; chomp $crispr_start_stop_line;
			$crispr_start_stop_line =~ s/\r//;
			$crispr_start_stop_line =~ s/^#//;
			$crispr_start_stop_line =~ s/ Crispr_begin_position: //;
			$crispr_start_stop_line =~ s/ Crispr_end_position: //;
								
			my @arr_tmp=split('\t', $crispr_start_stop_line);
			$rangeStart = $arr_tmp[0];
			$rangeEnd = $arr_tmp[1];
			$crispr_length = $rangeEnd - $rangeStart + 1;
								
			my $dr_line = $arr_spacer_file[$i + 8]; chomp $dr_line; 
			$dr_line =~ s/\r//; 
			$dr_line =~ s/# DR: //;
			
			my ($dr, $tmp1) = split('\t',$dr_line); 
			$reprep = $dr;
			
			my $crispr_id = "CRISPR$crispr_index" . "_" . $rangeStart . "_" . $rangeEnd;
			
			my $line0 = "$orgName\t$app\trepeat_region\t$rangeStart\t$rangeEnd\t$crispr_length\t$strand\t$xref\t";
			$line0 .= "ID=$crispr_id;Name=$crispr_id;Note=$reprep;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
			print WR "$line0\n";
			
			$i = $i + 11;
			my $repeatIndex = 0;
			my $spacerIndex = 0;
			
			my $x0 = $rangeStart;
			while($arr_spacer_file[$i] =~ / +([0-9]+)\t +([0-9]+)\t +([ACGT]+)/) {
				$spacerIndex++;
				$repeatIndex++;
				
				my $spacer_start = int($1);
				my $spacer_end = int($1) + int($2) - 1;
				my $spacer_length = int($2);
				my $spacer_seq = $3;
				my $spacer_id = "CRISPR$crispr_index" . "_" . "SPACER$spacerIndex" . "_" . $spacer_start . "_" . $spacer_end;
				
				my $repeat_start = $x0;
				my $repeat_end = $spacer_start - 1;
				my $repeat_length = $repeat_end - $repeat_start + 1;
				my $repeat_seq = $reprep;
				my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeatIndex" . "_" . $repeat_start . "_" . $repeat_end;
				
				my $line1 = "$orgName\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$repeat_length\t$strand\t$xref\t";
				$line1 .= "ID=$repeat_id;Name=$repeat_id;Parent=$crispr_id;Note=$repeat_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
				print WR "$line1\n";
				
				my $line2 = "$orgName\t$app\tbinding_site\t$spacer_start\t$spacer_end\t$spacer_length\t$strand\t$xref\t";
				$line2 .= "ID=$spacer_id;Name=$spacer_id;Parent=$crispr_id;Note=$spacer_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
				print WR "$line2\n";
				
				$x0 = $spacer_end + 1;
				
				if($arr_spacer_file[$i + 2] =~ /#====/) {
					$repeatIndex++;
					
					my $repeat_start = $x0;
					my $repeat_end = $rangeEnd;
					my $repeat_length = $repeat_end - $repeat_start + 1;
					my $repeat_seq = $reprep;
					my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeatIndex" . "_" . $repeat_start . "_" . $repeat_end;
				
					my $line1 = "$orgName\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$repeat_length\t$strand\t$xref\t";
					$line1 .= "ID=$repeat_id;Name=$repeat_id;Parent=$crispr_id;Note=$repeat_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=NA";
					print WR "$line1\n";
					last;
				}
				$i = $i + 2;	 
			}
		}
	}
	close(WR);
	return 1;
}

sub CCFinderToGFF {
	my ($path, $out) = @_;
	
	my $json;
	my $strand = "+";
	my $app = "CRISPRCasFinder";
	my $xref = ".";

	{
   		local $/;
    	open my $fh, '<', $path or die $!;
   		$json = <$fh>;
	}

	my $perl = decode_json $json;
	
	my @seqs = @{$perl->{Sequences}};
	
	open(WR,">$out");
	foreach my $seq (@seqs) {
		my $seq_id = $seq->{Id};	
		my @crisprs = @{$seq->{Crisprs}};
		
		my $count = 0;
		foreach my $crispr (@crisprs) {
			$count++;
			my $arr_start = $crispr->{Start};
			my $arr_end = $crispr->{End};
			my $score = $crispr->{Evidence_Level};
			my $ori = $crispr->{Potential_Orientation};
			my $reprep = $crispr -> {DR_Consensus};
			my $crispr_index = $count;
			my $crispr_id = "CRISPR$crispr_index" . "_" . $arr_start . "_" . $arr_end;
			my $crispr_length = $arr_end - $arr_start + 1;
			
			if ($ori eq "-") {
				$strand = "-";
			}
						
			my @regions = @{$crispr->{Regions}};
			
			my $spacer_count = 0;
			my $repeat_count = 0;
			
			my $line0 = "$seq_id\t$app\trepeat_region\t$arr_start\t$arr_end\t$crispr_length\t$strand\t$xref\t";
			$line0 .= "ID=$crispr_id;Name=$crispr_id;Note=$reprep;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=$score";
			print WR "$line0\n";
			
			foreach my $region (@regions) {
				my $type = $region->{Type};
				
				if ($type eq "DR") {
					$repeat_count++;
					my $repeat_start = $region->{Start};
					my $repeat_end = $region->{End};
					my $repeat_seq = $region->{Sequence};
					my $repeat_index = $repeat_count;
					my $repeat_length = length($repeat_seq);
					my $repeat_id = "CRISPR$crispr_index" . "_" . "REPEAT$repeat_index" . "_" . $repeat_start . "_" . $repeat_end;
					my $line1 = "$seq_id\t$app\tdirect_repeat\t$repeat_start\t$repeat_end\t$repeat_length\t$strand\t$xref\t";
				 	$line1 .= "ID=$repeat_id;Name=$repeat_id;Parent=$crispr_id;Note=$repeat_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=$score";
				 	print WR "$line1\n";
				}
				
				if ($type eq "Spacer") {
					$spacer_count++;
					my $spacer_start = $region->{Start};
					my $spacer_end = $region->{End};
					my $spacer_seq = $region->{Sequence};
					my $spacer_index = $spacer_count;
					my $spacer_length = length($spacer_seq);
					my $spacer_id = "CRISPR$crispr_index" . "_" . "SPACER$spacer_index" . "_" . $spacer_start . "_" . $spacer_end;
					my $line2 = "$seq_id\t$app\tbinding_site\t$spacer_start\t$spacer_end\t$spacer_length\t$strand\t$xref\t";
					$line2 .= "ID=$spacer_id;Name=$spacer_id;Parent=$crispr_id;Note=$spacer_seq;Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score=$score";
					print WR "$line2\n";
				}
			}
		}
	}
	close(WR);
	return 1;
}
