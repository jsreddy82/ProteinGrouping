#!/usr/bin/perl
#Program to group proteins identified with the same set of (identical) peptides
#Input file is Xtandem filtered output file from ProcessSearchResults.pl at : "http://sourceforge.net/projects/proteomicstools/"
use strict;
use Data::Dumper qw(Dumper);
if(@ARGV!=4)
{
	print "Program to group proteins identified with the same set of (identical) peptides in a filtered Xtandem output\n";
        print "\nUSAGE: \$perl group_Xtandem_Filtered.pl <Xtandem Filtered File> <Retained> <Protein-groups> <Processed Output>\n";
        print "\n<Xtandem Filtered File>: Output File from ProcessSearchResults.pl\n";
	print "<Retained>: Name of output file to save list of retained proteins in the processed Xtandem output\n";
	print "<Protein-groups>: Name of output file to save list of all proteins idenfied with Xtandem and their group IDs\n";
	print "<Processed Output>: Name of the output file to save the processed Xtandem Output retaining only one protein for each unique set of peptides (grouped)\n";
	exit 0;
}
open FILE, "$ARGV[0]" or die $!;
open OUT, "$ARGV[1]" or die $!;		#Final list of proteins that were retained in processed X!tandem_filtered file after grouping
open GOUT, "$ARGV[2]" or die $!; 	#List of proteins and their group ids.
open EXCL, "$ARGV[3]" or die $!; 	#Processed X!tandem_filtered output retaining unique list of identified proteins (Proteins identifed with identical peptides excluded, but redundant protein/s (ID) listed alongside retained protein from species of interest)  
my $line_count=1;
my $protein="";
my @lines=();
my %mapped=();
my @keep=();
my @data=();
my @grouped_all=();
my %prot_group=();
while(<FILE>)
{
	chomp($_);push(@data,$_);
	if($line_count<=12)
	{
		$line_count++;
		next;
	}
	if($_ =~ /^(Engine)/)
	{
		next;
	}
	elsif($_ =~ /^(Xtandem)/)
	{
		push (@lines,$_);
	}
	else
	{
		if(scalar(@lines)!=0)
		{
			my $i=0;
			foreach my $spectra(@lines)
			{
				my @peptide=split ("\t",$spectra);
				#$peptide[3]=~ s/(\+|\-)[0-9]*\.[0-9]*//g;
				$mapped{$protein}{$i}=join("\t",$peptide[1],$peptide[2],$peptide[3]);$i++;
			}
		}
		my @pid = split (/\t/,$_);
		my @pname = split (/\|/,$pid[0]);
		my @species= split (/\s/,$pname[2]);
		$protein= join("|",$pname[1],$species[0]);
		@lines=();	
	}
}close(FILE);
#print Dumper \%mapped;
my %temp_mapped=%mapped;
#print Dumper \%temp_mapped;
my $group=1;
my @included=();
foreach my $cprot(sort keys %mapped)
{
	if(grep {$_ eq $cprot} @included)
	{
		next;
	}
	
	my $matched=0;
	my @grouped=();
	my @temp=();
	push(@temp,$cprot);	
	foreach my $tprot(sort keys %temp_mapped)
	{
		if($cprot eq $tprot)
		{
			next;
		}
		elsif((keys %{$mapped{$cprot}}) == (keys %{$temp_mapped{$tprot}}))
		{
			if(grep {$_ eq $tprot} @included)
			{
				next;
			}
			else
			{
				my $count=0;	
				foreach my $cspectra(sort keys %{$mapped{$cprot}})
				{
					foreach my $tspectra(sort keys %{$temp_mapped{$tprot}})	
					{
						if($mapped{$cprot}{$cspectra} eq $temp_mapped{$tprot}{$tspectra})
						{
							#print "$mapped{$cprot}{$cspectra}--$temp_mapped{$tprot}{$tspectra}\n";
							$count++;
						}
					}
				}
				if($count == keys %{$temp_mapped{$tprot}})
				{
					push (@grouped,$tprot);
					$matched=1;
				}
			}
		}
	}
	if($matched==0)
	{
		push(@keep,$cprot);
	}
	my $flag=0;
	if(scalar @grouped > 0)
	{
		push (@grouped,$cprot);
		my @sorted_grouped=sort @grouped;
		foreach my $prot(@sorted_grouped)
		{
			$prot_group{$prot}=$group;
			print GOUT "$group\t$prot\n"; 
			if($prot =~ /GOSHI/)                                ############### CHANGE THIS KEYWORD TO SPECIES OF INTEREST
			{
				$flag=1;
				push (@keep,$prot);
			}
			push(@temp,$prot);push(@grouped_all,"$prot");
		}
		if($flag==0)
		{
			
			push(@keep,$sorted_grouped[0]);
		}			
		$group++;
	}
	foreach (@temp)
	{
		push (@included,$_);
	}	
} 
foreach my $id(@keep)
{
	print OUT "$id\n";				##### LIST OF RETAINED PROTEINS IN FILTERED FILE AFTER GROUPING
}
my $flag=0;
foreach my $line(@data)
{
	$line=~s/\r//g;
	if($line =~ /^(sp)|^(tr)/)
	{
		my @pid = split (/\t/,$line);
		my @pname = split (/\|/,$pid[0]);
		my @species= split (/\s/,$pname[2]);
		my $protein= join("|",$pname[1],$species[0]);
		#print $pname[1]."\n";
		if(grep {$_ eq $protein} @keep)
		{
			$flag=0; #print $protein."\n";
			if(defined $prot_group{$protein})
			{
				my $grp=$prot_group{$protein};
				my @gropd =();
				foreach my $id(sort keys %prot_group)
				{
					if ($prot_group{$id} ==$grp)
					{
						push(@gropd, $id);
					}
				}
				my $gids=join("\t",@gropd);
				$line=join("\t",$line,"GROUP$grp",$gids); ############### Add grouped protein ids to the current protein
			}
		}
		else
		{
			$flag=1;
		}
	}
	if($flag==0)
	{
		print EXCL $line."\n";
	}
}
close(OUT);
close(GOUT);
close(EXCL);
