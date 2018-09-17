=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Stephen Kazakoff <sh.kazakoff@gmail.com>
    
=cut

=head1 NAME

 PROVEAN

=head1 SYNOPSIS

 mv PROVEAN.pm ~/.vep/Plugins
 ./vep -i variants.vcf --plugin PROVEAN,/path/to/cache_dir,/path/to/jobs.txt

=head1 DESCRIPTION

 A VEP plugin to assist in the creation of the inputs needed to run
 PROVEAN. This plugin is also able to integrate PROVEAN scores if a cache
 directory has been populated and the PROVEAN output was saved using a
 '.out' extension.

 Basically, we run VEP to create a set of inputs for PROVEAN. We then run
 PROVEAN. Finally, we re-run VEP to incorporate the scores. We use an on-
 disk cache since it's not practical to pre-calculate all possible scores
 for insertions/deletions etc.


 1. Install the tools needed to run PROVEAN.

 1a. Install NCBI-BLAST:

 > wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-src.tar.gz
 > tar -xvf ncbi-blast-2.7.1+-src.tar.gz
 > cd ncbi-blast-2.7.1+-src/c++
 > ./configure --prefix=/software/ncbi-blast/ncbi-blast-2.7.1+
 > make
 > make install

 1b. Install CD-HIT:

 > wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz
 > mkdir -p /software/cd-hit/cd-hit-v4.6.8-2017-1208
 > tar -xvf cd-hit-v4.6.8-2017-1208-source.tar.gz -C /software/cd-hit
 > cd /software/cd-hit/cd-hit-v4.6.8-2017-1208
 > make

 1c. Install PROVEAN:

 > wget https://downloads.sourceforge.net/project/provean/provean-1.1.5.tar.gz
 > tar -xvf provean-1.1.5.tar.gz
 > cd provean-1.1.5
 > ./configure --prefix=/software/provean/provean-1.1.5 PSIBLAST=/ CDHIT=/ BLASTDBCMD=/
 > make
 > make install

 1d. Download the NCBI NR BLAST database:

 > wget --recursive --no-parent --accept 'nr.*' ftp://ftp.ncbi.nih.gov/blast/db
 > cd ftp.ncbi.nih.gov/blast/db
 > tar -xvf *.tar.gz


 2. Run VEP with the PROVEAN plugin for each of your samples. This will create the PROVEAN inputs:

 > mkdir -p /path/to/PROVEAN/cache_dir
 > vep -i variants.vcf --plugin PROVEAN,/path/to/PROVEAN/cache_dir,/path/to/PROVEAN/jobs.txt

 Note:
  * `cache_dir`: will contain subdirectories comprising a FASTA file and list of variants for each peptide.
  * `jobs.txt`: will contain a list of peptides that we will need to run PROVEAN for.


 3. Run the PROVEAN jobs and save the output into the cache.

 This example uses GNU Parallel on a single machine, but you may need to scale up if you have
 hundreds of thousands of peptides to analyse.

 3a. Install GNU Parallel if you don't already have it:

 > (wget -O - pi.dk/3 || curl pi.dk/3/) | bash

 3b. Run PROVEAN using GNU Parallel:

 > run_PROVEAN() {
 >   peptide="${1}"
 >   cache_dir="/path/to/PROVEAN/cache_dir"
 >
 >   query="${cache_dir}/${peptide}/${peptide}.fasta"
 >   variation="${cache_dir}/${peptide}/${peptide}.var"
 >   output="${cache_dir}/${peptide}/${peptide}.out"
 >
 >   provean.sh -q "${query}" -v "${variation}" > "${output}"
 > }
 > export -f run_PROVEAN
 > cat /path/to/PROVEAN/jobs.txt | parallel run_PROVEAN

 It is recommended that you also use the `--save_supporting_set` and `--supporting_set` options
 to skip subsequent BLAST searches / clustering steps. There's also a `--num_threads` option to
 reduce the BLAST search time.

 Please ensure all PROVEAN jobs have completed successfully.


 4. Re-run VEP to incorporate the scores from the cache directory:

 > vep -i variants.vcf --plugin PROVEAN,/path/to/PROVEAN/cache_dir,/path/to/PROVEAN/new_jobs.txt


 5. Repeat steps 3 and 4 if a file called `new_jobs.txt` was created.


=cut

package PROVEAN;

use strict;
use warnings;

use Fcntl qw(:flock :DEFAULT);

use Bio::SeqIO;
use Bio::SeqUtils;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    PROVEAN => "PROVEAN score"
  };
}

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $cache_dir = $self->params->[0];
  my $job_file = $self->params->[1];

  die("ERROR: PROVEAN cache dir not specified\n") unless $cache_dir;
  die("ERROR: PROVEAN job file not specified\n") unless $job_file;

  die("ERROR: PROVEAN cache dir not found\n") unless -d $cache_dir;
  die("ERROR: PROVEAN cache dir not writable\n") unless -w $cache_dir;

  $self->{cache_dir} = $cache_dir;
  $self->{job_file} = $job_file;

  return $self;
}

sub run {
  my ($self, $tva) = @_;

  my $hgvs_p = $tva->hgvs_protein;

  return {} unless $hgvs_p;

  my ($protein_id, $variation) = $hgvs_p =~ /^(.*?):p\.([^:]+)$/;

  return {} unless $protein_id && $variation;

  my $ref_seq = $tva->transcript_variation->_peptide;

  # ignore silent changes
  return { PROVEAN => 0 } if $variation =~ /=$/;

  my %onecode = %Bio::SeqUtils::ONECODE;

  # translate all three letter codes to one letter codes
  $variation =~ s/$_/$onecode{$_}/g for keys %onecode;

  # ensure simple deletions are handled correctly
  $variation =~ s/\*$/del/ unless $variation =~ /[a-z]/;

  my $dir = join('/', $self->{cache_dir}, $protein_id);

  mkdir $dir unless -d $dir;

  my $fasta_file = "$dir/$protein_id.fasta";
  my $var_file = "$dir/$protein_id.var";
  my $output_file = "$dir/$protein_id.out";

  # create a FASTQ file if it doesn't exist
  if (sysopen(my $fasta_fh, $fasta_file, O_CREAT | O_EXCL | O_WRONLY)) {

    my $seqIO = Bio::SeqIO->new(-fh => $fasta_fh, -format => "fasta");

    my $seq = Bio::Seq->new(-seq => $ref_seq, -id => $protein_id);

    $seqIO->write_seq($seq);

    close $fasta_fh;
  }

  my %provean_scores;

  # try to read any PROVEAN results 
  if (sysopen(my $output_fh, $output_file, O_RDONLY)) {

    while (<$output_fh>) {
      chomp;

      next if /^#/ || !/\t/;

      my ($provean_variation, $provean_score) = split /\t/;

      $provean_scores{$provean_variation} = $provean_score;
    }

    close $output_fh;
  }

  my $f;

  # try to add any novel variations to the variations file
  if (sysopen(my $var_fh, $var_file, O_CREAT | O_RDWR)) {

    flock($var_fh, LOCK_EX);

    my %vars = map { chomp; $_ => 1 } <$var_fh>;

    my %h;

    for ($variation, keys %provean_scores) {

      $h{$_} = 1 unless exists $vars{$_};
    }

    if (%h) {

      %vars = (%vars, %h);

      seek($var_fh, 0, 0);

      truncate($var_fh, 0);

      print $var_fh "$_\n" for sort keys %vars;

      $f = 1;
    }

    close $var_fh;
  }

  # if we added new variations, add a job to our job list
  if ($f) {

    my $job_file = $self->{job_file};

    if (sysopen(my $job_fh, $job_file, O_CREAT | O_RDWR)) {

      flock($job_fh, LOCK_EX);

      my %jobs = map { chomp; $_ => 1 } <$job_fh>;

      print $job_fh "$protein_id\n" unless exists $jobs{$protein_id};

      close $job_fh;
    }
  }

  my $score = $provean_scores{$variation};

  return $score ? { PROVEAN => $score } : {};
}

1;

