=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

  Ensembl <https://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

  PrimateAI

=head1 SYNOPSIS

  mv PrimateAI.pm ~/.vep/Plugins

  ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz
  ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz

=head1 DESCRIPTION

  The PrimateAI VEP plugin is designed to retrieve clinical impact scores of variants, as described in https://www.nature.com/articles/s41588-018-0167-z. Please consider citing the paper if using this plugin.
  In brief, common missense mutations in non-human primate species are usually benign in humans. Thousands of common variants from six non-human primate species were used to train a deep neural network to identify pathogenic mutations in humans with a rare disease.

  This plugin uses files generated by the PrimateAI software, which is available from https://github.com/Illumina/PrimateAI. The files containing predicted pathogenicity scores can be downloaded from https://basespace.illumina.com/s/yYGFdGih1rXL (a free BaseSpace account may be required):
      PrimateAI_scores_v0.2.tsv.gz (for GRCh37/hg19)
      PrimateAI_scores_v0.2_hg38.tsv.gz (for GRCh38/hg38)

  Before running the plugin for the first time, the following steps must be taken to format the downloaded files:
  1.  Unzip the score files
  2.  Add '#' in front of the column description line
  3.  Remove any empty lines.
  4.  Sort the file by chromosome and position
  5.  Compress the file in .bgz format
  6.  Create tabix index (requires tabix to be installed).

  Command line examples for formatting input files:
    > gunzip -cf PrimateAI_scores_v0.2.tsv.gz | sed '12s/.*/#&/' | sed '/^$/d' | awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | bgzip > PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz
    > tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz

    > gunzip -cf PrimateAI_scores_v0.2_hg38.tsv.gz | sed '12s/.*/#&/' | sed '/^$/d' | awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | bgzip > PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz
    > tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz

=cut

package PrimateAI;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $file = $self->params->[0];
  die("ERROR: PrimateAI scores file $file not found\n") unless $file && -e $file;

  my $assembly = $self->{config}->{assembly};

  $self->add_file($file);

  #Check that the file matches the assembly
  if ($file =~ /GRCh37/){
    if ($assembly ne "GRCh37"){
      die "The PrimateAI scores file contains GRCh37 coordinates, but you have not selected the GRCh37 assembly in your input. Please check and try again.\n";
    }
  } elsif ($file =~ /GRCh38/){
      if ($assembly ne "GRCh38"){
        die "The PrimateAI scores file contains GRCh38 coordinates, but you have not selected the GRCh38 assembly in your input. Please check and try again.\n";
      }
  } else{
      print "Using $file with assembly $assembly\n";
  }

  $self->expand_left(0);
  $self->expand_right(0);

  return $self;
}

sub feature_types {
  return ['Feature', 'Intergenic'];
}

sub get_header_info {
  my $self = shift;

  return {
    PrimateAI => "PrimateAI score for variants"
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $allele = $tva->variation_feature_seq;

  return {} unless $allele =~ /^[ACGT]$/;

  #Get the start and end coordinates, and ensure they are the right way round (i.e. start < end).
  my $vf_start = $vf->{start};
  my $vf_end = $vf->{end};
  ($vf_start, $vf_end) = ($vf_end, $vf_start) if $vf_start > $vf_end;

  #Check the strands and complement the allele if necessary.
  if ($vf->{strand} <0){
    reverse_comp(\$allele);
  }

  #Compare the position and allele
    my ($res) = grep {
    $_->{start} == $vf_start &&
    $_->{end} == $vf_end &&
    $_->{alt} eq $allele
  } @{$self->get_data($vf->{chr}, $vf_start, $vf_end)};

  #Return data if matched
  return $res ? $res->{result} : {};
}

sub parse_data {
  my ($self, $line) = @_;

  #Necessary columns from the input file.
  my @line = split(/\t/, $line);
  my ($chr, $pos, $alt, $score) = ($line[0], $line[1], $line[3], $line[10]);

  return {
    alt => $alt,
    start => $pos,
    end => $pos,
    result => {
      #Score to be returned
      PrimateAI   => $score
    }
  };
}

sub get_start {
  return $_[1]->{'start'};
}

sub get_end {
  return $_[1]->{'end'};
}

1;
