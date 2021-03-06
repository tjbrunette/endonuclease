#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

# script for running psi-blast and keeping the output around.

my $blastpgp    = "/work/tex/src/blast/bin/blastpgp";
my $e_value     = 1e-15;
my $h_value     = 1e-15;
my $n_rounds    = 1;
my $n_procs     = 20;
my $archive     = 0;
my $database    = '/scratch/ROSETTA/genomes/nr';
#my $database     = '/work/brunette/databases/endo/mito.aa';

if ( exists $ENV{PSIBLAST_NPROCS} ) {
	$n_procs = $ENV{PSIBLAST_NPROCS};
}

GetOptions(
	"e_value=f"  => \$e_value,
	"h_value=f"  => \$h_value,
	"n_rounds=i" => \$n_rounds,
	"n_procs=i"  => \$n_procs,
	"archive"    => \$archive,
	"db=s"       => \$database,
	"blastpgp=s" => \$blastpgp,
);

my $usage = <<USAGE;
usage: $0 [options] fasta1 fasta2

options are:
	--e_value  <e_value>  (threshold for inclusion of sequenced in output alignments, defaults to 0.001)
	--h_value  <h_value>  (threshold for inclusion of sequences in profile, defaults to 0.001)
	--n_rounds <n_rounds> (run this many rounds of psi-blast, defaults to 10)
	--n_procs  <n_procs>  (use this many processors)
	--archive             (put all of the files produced into fasta.tar.gz)
	--db       <db>       (path to blast database)
	--blastpgp <blastpgp> (path to blastpgp binary)
USAGE

my @fasta_files = @ARGV;
if ( scalar(@fasta_files) == 0 ) {
	warn $usage;
	exit 1;
}

foreach my $fasta (@fasta_files) {
	if ( ! -f $fasta ) { die "Error: can't open file $fasta!\n$usage\n"; }

	my @files_created;

	# psiblast round 1
	my $round = 1;
	my $chk_file  = "$fasta.$round.chk";
	my $aln_file  = "$fasta.$round.psiblast";
	my $pssm_file = "$fasta.$round.pssm";
	my $cmd = "$blastpgp -i $fasta -d $database -o $aln_file -a $n_procs -C $chk_file -j2 -Q $pssm_file -e $e_value -h $h_value";
	#print $cmd, "\n";
	system( $cmd );
	push @files_created, $chk_file;
	push @files_created, $aln_file;
	push @files_created, $pssm_file;

	for my $round ( 2 .. $n_rounds ) {
		my $last_round    = $round - 1;
		my $aln_file      = "$fasta.$round.psiblast";
		my $chk_file      = "$fasta.$round.chk";
		my $last_chk_file = "$fasta.$last_round.chk";
		my $pssm_file     = "$fasta.$round.pssm";

		my $cmd = "$blastpgp -i $fasta -d $database -R $last_chk_file -o $aln_file -a $n_procs -C $chk_file -Q $pssm_file -j2";
		#print $cmd, "\n";
		system( $cmd );

		push @files_created, $chk_file;
		push @files_created, $aln_file;
		push @files_created, $pssm_file;
	}

	if ( $archive ) {
		my $tar_file   = "$fasta.tar.gz";
		my $dir        = dirname  $tar_file;
		$tar_file      = basename $tar_file;
		my $rm_cmd     = "rm " . join ' ', @files_created;

		@files_created = map { basename $_ } @files_created;
		my $tar_cmd    = "cd $dir; tar czvf $tar_file " . join ' ', @files_created;

		if ( -f $tar_file ) {
			warn "Error: not overwriting file $tar_file.\n";
			next;
		}

		print $tar_cmd, "\n";
		system( $tar_cmd );
		system( $rm_cmd );
	}
}
