#!/usr/bin/perl

use strict;
use warnings;

#use lib '/work/tex/src/cm_scripts/perl_lib';
use FindBin;
use lib "$FindBin::Bin/../perl_lib";

use Cwd qw/ getcwd abs_path /;
use Getopt::Long;
use File::Basename;
use File::Copy;
use Rosetta::Job;
use Rosetta::Util qw/ options_to_str trim_whitespace script_dir /;
use hhsearch::Util;

my %options;
$options{script_base}   = script_dir();
$options{n_procs}       = 1;
$options{hhsearch_path} = '/work/tex/hhsearch';
$options{n_rounds}      = 3;
$options{use_hhblits}   = 0;
$options{talos_fn}      = 0;

&GetOptions(
	\%options,
	"db=s",
	"n_procs=i",
	"n_rounds=i",
	"use_hhblits!",
	"hhsearch_path=s",
	"talos_fn=s",
);

print STDERR "running with the following options:\n";
print STDERR options_to_str( \%options ), "\n";

my @fns = @ARGV;

foreach my $fn (@fns) {
	my $hhm_file;

	if ( ! -f $fn ) {
		warn "Error: file $fn doesn't exist!\n";
		next;
	}

	if ( abs_path(getcwd) ne abs_path(dirname($fn)) ) {
		copy $fn, getcwd;
	}
	$fn = basename $fn;
	if ( $fn =~ /\.fasta$/ ) {
		$hhm_file = make_hhm( $fn, '.', \%options );
	} elsif ( $fn =~ /\.pdb$/ ) {
		system( "$options{script_base}/pdb2fasta.py $fn > $fn.fasta" );
		$hhm_file = make_hhm( "$fn.fasta", '.', \%options );
		#add_dssp_to_hhm($fn,$hhm_file,\%options);
	}
}

sub add_dssp_to_hhm {
	my $pdb_file = shift;
	my $hhm_file = shift;
	my $options  = shift;

	my $dssp_fn = "$pdb_file.ss_dssp";
	if ( ! -f $dssp_fn ) {
		system( "$options->{script_base}/dssp.pl $pdb_file" );
	}

	open FILE, "<$dssp_fn" or die $!;
	my $dssp_ss = '';
	while ( my $line = <FILE> ) {
		if ( $line !~ /^>/ ) {
			chomp $line;
			$dssp_ss .= trim_whitespace($line);
		}
	}
	close FILE or die $!;

	$dssp_ss =~ s/L/C/g;

	open FILE, "<$hhm_file" or die $!;
	my @file = <FILE>;
	close FILE or die $!;

	open FILE, ">$hhm_file" or die $!;
	foreach my $line (@file) {
		if ( $line =~ /^SEQ/ ) {
			print FILE $line;
			print FILE ">ss_dssp\n$dssp_ss\n";
		} else {
			print FILE $line;
		}
	}
	close FILE or die $!;
}
