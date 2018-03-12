#!/usr/bin/perl

use strict;
use warnings;

#script to report ATAC-seq summary from outputs after running Ting's pipeline "run_ATACseq.sh"
#Ting Wang
#2016.07

if(@ARGV != 2){
	die "wrong inputs\nUSAGE :: perl summarize_RNAseq.pl <directory> <outfile>\n";
}else{
	open my $outfile, ">$ARGV[1]" or die "wrong outfile\n";

	print $outfile "Sample\tTotal_Read_pairs\tTrimmed_Read_pairs\tRemaining_Read_pairs\tUniquely-concordantly-mapped_Rate\tMulti-concordantly-mapped_Rate\tOverall-mapped_Rate\tchrM.LowMAPQ_Rate\tDuplicate_Rate\tBlacklist_Rate\tRemaining_Reads\tPeaks\n";
	
	my $dir = $ARGV[0];

	#do rest of stuff
	opendir(DIR, $dir) or die $!;
	foreach my $file (sort readdir(DIR)){
		next unless (-d "$dir/$file");

		#tried "ne" and it didn't work (?)
		if($file eq "." || $file eq ".."){}else{
			#1 name
			print $outfile "$file\t";

			#total reads and trimmed reads and remaining
			my $tot_reads = 0;
                        my $trimmed_reads = 0;
			my $remain_reads = 0;
			opendir(SUBDIR, "$dir/$file");
			my @trim_report_file = grep(/trimming_report.txt$/, sort readdir(SUBDIR));
			open my $trim_report, "$dir/$file/$trim_report_file[1]" || die "trimming_report.txt open error\n";
			my $numlines = () = <$trim_report>;
			seek($trim_report, 0, 0);
			
			#my $tot_reads = 0;
			#my $trimmed_reads = 0;

			my $lc = 1;
			while(<$trim_report>){
				if($lc == $numlines - 4){
					my @line = split /\s+/;
					$tot_reads = $line[0];
				}
				if($lc == $numlines ){
					my @line = split /\s+/;
					$trimmed_reads = $line[-2];
				}	
				$lc++;
			}

			$remain_reads = $tot_reads - $trimmed_reads;						
			print $outfile "$tot_reads\t$trimmed_reads\t$remain_reads\t";	
			close($trim_report);
		
			#mapping rates
			open my $log_final, "$dir/$file/bowtie2align.log" or die "bowtie2align.log error\n";
			my $line8 = 0;
			my $line9 = 0;
			while (<$log_final>){
				if($. == 8){
					my @line = split /\s+/;
					$line8 = $line[2];
					$line8 =~ s/\(//g;
					$line8 =~ s/\)//g;
					print $outfile "$line8\t";
				}

				if($. == 9){
					my @line = split /\s+/;
					$line9 = $line[2];
					$line9 =~ s/\(//g;
					$line9 =~ s/\)//g;
					print $outfile "$line9\t";
				}

				if($. == 19){
					my @line = split /\s+/;
					print $outfile "$line[0]\t";	
				}
			}
			close($log_final);

			#chrM rate
			my $read_rm_chrM = 0;
			my $chrM_rate = 0;
			open my $rmchrM_file, "$dir/$file/bowtie2align.nochrM.sorted.flagstat" or die "rm chrM file error\n";
			while(<$rmchrM_file>){
				chomp;
				my @line = split /\s+/;
				if($. == 1){
					$read_rm_chrM = $line[0];
				}
			}
			$chrM_rate = ($remain_reads * 2 - $read_rm_chrM) / $remain_reads / 2 * 100;
			my $rounded = sprintf("%.2f",$chrM_rate);
			print $outfile "$rounded"."%\t";
			close($rmchrM_file);

			#duplicate rate
			my $read_rmdup = 0;
			my $dup_rate = 0;
			open my $rmdup_file, "$dir/$file/bowtie2align.nochrM.sorted.rmdup.flagstat" or die "rmdup file error\n";
                        while(<$rmdup_file>){
                                chomp;
                                my @line = split /\s+/;
                                if($. == 1){
                                        $read_rmdup = $line[0];
                                }
                        }
                        $dup_rate = ($read_rm_chrM - $read_rmdup) / $read_rm_chrM * 100;
                        my $rounded2 = sprintf("%.2f",$dup_rate);
                        print $outfile "$rounded2"."%\t";
			close($rmdup_file);

			#blacklist rate
			my $read_rmbl = 0;
			my $bl_rate = 0;
			open my $rmbl_file, "$dir/$file/bowtie2align.nochrM.sorted.rmdup.rmblacklist.flagstat" or die "rm blacklist file error\n";
                        while(<$rmbl_file>){
                                chomp;
                                my @line = split /\s+/;
                                if($. == 1){
                                        $read_rmbl = $line[0];
                                }
                        }
                        $bl_rate = ($read_rmdup - $read_rmbl) / $read_rmdup * 100;
                        my $rounded3 = sprintf("%.2f",$bl_rate);
                        print $outfile "$rounded3"."%\t";
			close($rmbl_file);

			#remaining mates
			print $outfile "$read_rmbl\t";

			#peaks
			my $peaks = 0;
			open my $peak_file, "$dir/$file/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_peaks.narrowPeak" or die "pak file error\n";
			while(<$peak_file>){
				$peaks++;
			}
			print $outfile "$peaks";

			#all done
			print $outfile "\n";
		}
	}

	closedir(DIR);
}

#transform matrix
open(FD, "$ARGV[1]") || die;
my $line = <FD>;
my @out = split(/\s+/, $line);
while($line = <FD>){
	my @arr = split(/\s+/, $line);
	for (my $i = 0; $i < @arr; $i++){
		$out[$i].="\t".$arr[$i];
	}
}
close(FD);
open(FD,">$ARGV[1].tmp") || die;
for (my $i = 0; $i < @out; $i++){
	print FD $out[$i]."\n";
}
close(FD);
system("mv $ARGV[1].tmp $ARGV[1]");
exit 0;

