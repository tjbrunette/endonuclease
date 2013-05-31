package Rosetta::Template;

use Exporter;
@ISA = qw/ Exporter /;

@EXPORT_OK = qw/
	ids_from_aln_file
	filter_alignments_by_id
	make_vall_lines
/;

use Cwd qw/ abs_path /;
use Rosetta::PDB;
use Rosetta::Job;
use Rosetta::Util;
use Rosetta::JobList;
use Rosetta::Sequence;

use blast;
use alignment;

use List::Util qw/ min /;
use File::Basename;

sub ids_from_aln_file {
	my $aln_fn = shift;
	my @alns = @{ alignment::parse_alignments( $aln_fn ) };

	my @template_ids = map { substr( $_, 0, 5 ) } map { $_->template_name }
		@alns;

	# make a unique list of template_ids
	my %count;
	my @unique;
	foreach my $id (@template_ids) {
		if ( !exists $count{$id} ) {
			$count{$id}++;
			push @unique, $id;
		}
	}

	return \@unique;
	#return unique( \@template_ids );
}

sub make_vall_lines {
	my $pdb_file   = shift;
	my $options    = shift;

	my $dir        = dirname($pdb_file);
	$pdb_file      = basename($pdb_file);
	my $fasta_file = "$pdb_file.fasta";

	if ( ! -f "$dir/$pdb_file" ) {
		print STDERR "No PDB file found ($dir/$pdb_file).\n";
		exit 1;
	}
	system( "pdb2fasta.py $dir/$pdb_file > $dir/$fasta_file" );

	if ( ! -f "$dir/$fasta_file" ) {
		print STDERR "No fasta file found ($dir/$fasta_file).\n";
		exit 1;
	}

	my $ideal_dir = assemble_path( $dir, 'ideal' );
	my $vall_dir  = assemble_path( $dir, 'vall' );
	my $pssm_dir  = assemble_path( $dir, 'pssm' );

	my $pssm_job = Rosetta::Sequence::pssm_from_fasta(
		assemble_path($dir,$fasta_file), $pssm_dir,
	);
	my $idealize_job = Rosetta::PDB::idealize_pdb(
		assemble_path($dir,$pdb_file), $ideal_dir, $options
	);
	$pssm_job->run_with_message( "making pssm for $fasta_file" );
	$idealize_job->run_with_message( "making ideal pdb for $pdb_file" );

	# extract torsions
	my $instance = Rosetta::Instance->new(
		mini_prefix => $options->{mini_path},
		mode        => $options->{mini_compile_mode},
		compiler    => $options->{mini_compiler},
		db_path     => $options->{mini_db_path},
	);
	my $results_fn = $idealize_job->results_file;
	my $job = $instance->generate_job( "angles", "-mute all -in:file:s $results_fn > $results_fn.torsions" );
	$job->results_file( "$results_fn.torsions" );
	$job->lockfile( "$results_fn.torsions.lock" );
	$job->dir( $ideal_dir );
	$job->run_with_message( "reading torsions from $results_fn" );

	my $vall_lines;

	my %torsions;
	open ANGLES, "<$ideal_dir/$results_fn.torsions" or die $!;
	<ANGLES>;
	while ( my $line = <ANGLES> ) {
		chomp $line;
		$line =~ s/^\s+//g;
		if ( $line =~ /^\d+\s+\w\s+\w\s+[\-\d\.]+/ ) {
			my ($idx,$res,$ss,$phi,$psi,$omega,$CA_x,$CA_y,$CA_z,$chi1) = split /\s+/, $line;
			$torsions{$idx} = {
				aa    => $res,
				resi  => $resi,
				ss    => $ss,
				phi   => $phi,
				psi   => $psi,
				omega => $omega,
				CA_x  => $CA_x,
				CA_y  => $CA_y,
				CA_z  => $CA_z,
				chi1  => $chi1,
			};
		} else {
			#warn "didn't match line $line\n";
		}
	}
	close ANGLES or die $!;

	my %profile;
	open PROFILE, "<$pssm_dir/$fasta_file.3.pssm" or die $!;
	print "reading $pssm_dir/$fasta_file.3.pssm\n";
	while ( my $line = <PROFILE> ) {
		chomp $line;
		$line =~ s/^\s+//g;
		if ( $line =~ /^\d+\s+\w\s+[\-\d]+/ ) {
			my ($res,$aa,@profile) = split /\s+/, $line;
			@profile = @profile[0 .. 19];
			$profile{$res} = \@profile;
		}
	}
	close PROFILE or die $!;

	my $pdb_id = lc(substr(basename($pdb_file),0,4)) . uc(substr(basename($pdb_file),4,1));

	my $vall_lines = '';
	use List::Util qw/ max /;
	my $max_res = max( keys %torsions );
	foreach my $resi ( sort { $a <=> $b } keys %torsions ) {
		my $nalign = 10;
		$vall_lines .= sprintf (
			"%5s %1s %1s %5d %4d %4d " .
			"%8.2f %8.2f %8.2f " .       # xyz
			"%8.3f %8.3f %8.3f %8.3f " . # torsions
			"%3d %4.2f %5.3f ",
			$pdb_id,
			$torsions{$resi}->{aa},
			$torsions{$resi}->{ss},
			$resi,
			0,
			0,
			$torsions{$resi}->{CA_x},
			$torsions{$resi}->{CA_y},
			$torsions{$resi}->{CA_z},
			$torsions{$resi}->{phi},
			$torsions{$resi}->{psi},
			$torsions{$resi}->{omega},
			$torsions{$resi}->{chi1},
			$nalign,
			$acc,
			$gap
		);

		my $new_prof = boltzmann_weight( $profile{$resi} );

		$vall_lines .= join ' ', ( @$new_prof );
		$vall_lines .= "\n";
	}
	open FILE, ">$ideal_dir/$pdb_id.vall" or die $!;
	print FILE $vall_lines;
	close FILE or die $!;
}

sub add_frags_from_templates {
	my $frag_fn = shift;
	my $aln_fn  = shift;
}

sub boltzmann_weight {
	my $p = shift;

	my $temp = 0.5;
	my $sum  = 0;
	foreach my $n (@$p) {
		$sum += exp( $temp * $n );
	}

	my @new_prof;
	foreach my $n (@$p) {
		push @new_prof, sprintf( "%0.3f", (1 / $sum * exp( $temp * $n )) );
	}

	return \@new_prof;
}

sub make_template_alignments {
	my $fasta_fn = shift;
	my $aln_fn   = shift;
	my $dir      = shift;
	my $options  = shift;

	my $aln_maker = $options->{make_alignments};
	my $options_dump = '';
	if ( exists $options->{use_hhblits} ){
		$options_dump = "-use_hhblits";
	}
	if( exists $options->{fast_alignment_gen} ){
		$options_dump = $options_dump . " -fast_alignment_gen";
	}
	mkdir_safe( $dir );
	copy_safe( $fasta_fn, $dir );
	my $job = Rosetta::Job->new(
		executable   => $aln_maker,
		args         => [
			basename($fasta_fn),
			'-outfile', $aln_fn,
			'-min_alns', 50, $options_dump,
			'-max_template_pct_id', $options->{max_template_pct_id},
			'-aln_dir', $options->{aln_dir},
			'-targetid', $options->{target_id},
		],
		lockfile     => "make_alns.lock",
		dir          => $dir,
		results_file => $aln_fn,
	);
	return $job;
}

sub setup_templates {
	my $template_dir    = shift;
	my $aln_file        = shift;
	my $options         = shift;
	my $native_fn       = shift;

	my $new_aln_fn = join '.', ( basename($aln_file), "valid" );
	if ( -f assemble_path( $options->{aln_dir}, $new_aln_fn ) ) {
		return;
	}

	$options->{logger}( "setting up templates in directory $template_dir" );

	mkdir_safe( $template_dir );

	my @valid_ids;
	my @template_ids = @{ ids_from_aln_file( $aln_file ) };

	$options->{logger}( "validating templates!" );
	ID: foreach my $id (@template_ids) {
		my $pdbid = lc( substr( $id, 0, 4 ) );
		my $chain = uc( substr( $id, 4, 1 ) );
		my $pdb_file = "$pdbid$chain.pdb";
		$options->{logger}( "looking for $pdb_file ... " );
		if ( -f assemble_path( $template_dir, $pdb_file ) ) {
			#print " already exists as ", assemble_path($template_dir,$pdb_file), "\n";
			$options->{logger}(
				join ' ', ( " already exists as", assemble_path($template_dir,$pdb_file) )
			);
		} else {
			$options->{logger}( "attempting to fetch using get_pdb.py ... " );
			my $get_pdb = assemble_path( script_dir(), 'get_pdb.py' );

			my $cmd = "cd $template_dir; $get_pdb $pdbid $chain";
			my $output = `$cmd`;
			$options->{logger}($output);
			if ( ! -f assemble_path( $template_dir, $pdb_file ) ) {
				#print "Error getting pdb: $pdbid $chain\n";
				#print "output:\n$output\n";
				$options->{logger}( "Error getting pdb: $pdbid $chain" );
				$options->{logger}( "output:\n$output" );
			}
		}

		if ( ! -f assemble_path( $template_dir, $pdb_file ) ) {
			$options->{logger}( "error getting pdb: $pdbid $chain!" );
			next ID;
		} else {
			$options->{logger}( " found pdb for $pdbid $chain." );
		}

		my $valid_id = "$pdbid$chain";
		push @valid_ids, $valid_id;

		# make vall lines for each template
		#make_vall_lines( abs_path("$template_dir/$pdbid$chain.pdb"), $options );
	} # foreach my $id (@template_ids)

	# only include alignments that have valid pdb files
	$options->{logger}( "filtering alignments with these ids:" );
	$options->{logger}( @valid_ids );
	filter_alignments_by_id(
		basename($aln_file), $new_aln_fn, dirname($aln_file), \@valid_ids
	);
	my $aln_dir = dirname($aln_file);

	unlink(assemble_path($aln_dir,'alignment.filt'));
	make_local_link( $aln_dir, $new_aln_fn, 'alignment.filt' );
}


sub filter_alignments_by_id {
	my $aln_fn = shift;
	my $filtered_aln_file = shift;
	my $dir = shift;
	my $ids_to_include = shift;

	if ( -f assemble_path( $dir, $filtered_aln_file ) ) {
		print "Skipping run because $filtered_aln_file exists!\n";
		return;
	}

	my $start_dir = getcwd;
	chdir( $dir );
	my @alns = @{ alignment::parse_alignments( $aln_fn ) };

	my %id_count;
	foreach my $id (@$ids_to_include) {
		$id_count{$id} = 0;
	}

	my $n_printed = 0;
	open FILE, ">$filtered_aln_file" or die
		"Error opening file $filtered_aln_file ($!)";
	foreach my $aln (@alns) {
		$aln->template_name(
			lc( substr( $aln->get_template_name, 0, 4 ) ) .
			uc( substr( $aln->get_template_name, 4, 1 ) ) .
			  ( substr( $aln->get_template_name, 5 )    )
		);
		my $id = substr( $aln->template_name, 0, 5 );

		if ( exists $id_count{$id} ) {
			print FILE $aln->filt_string;
			$id_count{$id}++;
		}
	}
	close FILE or die $!;

	foreach my $id ( sort keys %id_count ) {
		my $count = $id_count{$id};
		print join ' ', ( "included $count alignments to", $id );
		print "\n";
	}

	chdir( $start_dir );
}

1;
