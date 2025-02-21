=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

 SpliceAI

=head1 SYNOPSIS

 mv SpliceAI.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_snv_.vcf.gz,
 indel=/path/to/spliceai_indel_.vcf.gz

=head1 DESCRIPTION

 A VEP plugin that retrieves pre-calculated annotations from SpliceAI.
 SpliceAI is a deep neural network, developed by Illumina, Inc 
 that predicts splice junctions from an arbitrary pre-mRNA transcript sequence.
 By default, this plugin appends all scores from SpliceAI files.

 Delta score of a variant, defined as the maximum of (DS_AG, DS_AL, DS_DG, DS_DL), 
 ranges from 0 to 1 and can be interpreted as the probability of the variant being 
 splice-altering. The author-suggested cutoffs are:
   * 0.2 (high recall)
   * 0.5 (recommended)
   * 0.8 (high precision)

 This plugin is available for both GRCh37 and GRCh38.

 More information can be found at:
 https://pypi.org/project/spliceai/

 Please cite the SpliceAI publication alongside VEP if you use this resource:
 https://www.ncbi.nlm.nih.gov/pubmed/30661751

 Running options:

  cutoff       : Only return the scores for the speficied cutoff
                 Accepted values are between 0 and 1

  split_output : Return each type of score in a different header.
                 This is easier for parsing the output file.

 Output: 
  The output includes the gene symbol, delta scores (DS) and delta positions (DP)
  for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL).

  - For tab the output contains one header "SpliceAI_pred" with all
    the delta scores and positions. The format is:
      "SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"

  - For JSON the output is a hash with the following format:
    "spliceai":
      {"DP_DL":0,"DS_AL":0,"DP_AG":0,"DS_DL":0,"SYMBOL":"X","DS_AG":0,"DP_AL":0,"DP_DG":0,"DS_DG":0}

  - For VCF output and option 'split_output' the delta scores and positions are stored in different headers.
    The values are "SpliceAI_pred_xx" being "xx" the score/position.
      Example: "SpliceAI_pred_DS_AG" is the delta score for acceptor gain.

  Gene matching:
  SpliceAI can contain scores for multiple genes that overlap a variant,
  and VEP can also predict consequences on multiple genes for a given variant.
  The plugin only returns SpliceAI scores for the gene symbols that match (if any).

 If plugin is run with option 2, the output also contains a flag: "PASS" if delta score
 passes the cutoff, "FAIL" otherwise. 

 The following steps are necessary to run this plugin:

 1. Download the input files:

  - The Illumina-generated files with the annotations for all possible substitutions (snv), 1 base insertions 
  and 1-4 base deletions (indel), within genes are available through basespace 
  (https://basespace.illumina.com/s/otSPW8hnhaZR). 
    - To download via Illumina's basespace:
      1. Log-in to your Illumina account or sign-up (for free) if you do not have one. 
      2. Once you're in, a "Share Project" pop-up will appear - click "accept". 
      3. A smaller pop-up in the bottom right will then read "Share Accepted". Click "Predicting splicing from primary sequence".
      4. You will get a list of files. Select "genome_scores_v1.3".
      5. You will get an info/landing page. Under "Analysis: genome_scores_v1.3", select "FILES". 
      6. Click the file icon next to "genome_scores_v1.3" and you will get a list of available files. 
      7. Click filenames to download the relevant files - note that raw/masked, hg19/hg38 and snv/indel files are available. 
  
  - The Ensembl-generated files with the annotations for all possible substitutions (snv), 1 base insertions, within genes are available through Ensembl 
 (https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/). Ensembl does not provide indel annotations, however, 
 Ensembl-generated files include annotations for Ensembl MANE select transcripts for v107 and v110 releases. 

 2. Tabix the files (if derived from Illumina). .tbi files are provided for Ensembl-derived VCFs. 

  - GRCh37:
  tabix -p vcf spliceai_scores.raw.snv.hg19.vcf.gz
  tabix -p vcf spliceai_scores.raw.indel.hg19.vcf.gz

  - GRCh38:
  tabix -p vcf spliceai_scores.raw.snv.hg38.vcf.gz
  tabix -p vcf spliceai_scores.raw.indel.hg38.vcf.gz

 3. The plugin can then be run: 

 - With Illumina files: 
  ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz
  ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5
  ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz,split_output=1

  - Or with Ensembl files: 
  ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_scores.masked.snv.ensembl_mane.grch38.110.vcf.gz,indel=/path/to/spliceai_scores.masked.indel.hg38.vcf.gz

=cut

package SpliceAI;

use strict;
use warnings;
use List::Util qw(max);

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;
use Bio::EnsEMBL::Variation::VariationFeature;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my $output_vcf;

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  $self->add_file($param_hash->{snv});
  $self->add_file($param_hash->{indel});

  if(defined($param_hash->{cutoff})) {
    my $cutoff = $param_hash->{cutoff};
    if($cutoff < 0 || $cutoff > 1) {
      die("ERROR: Cutoff score must be between 0 and 1!\n");
    }
    $self->{cutoff} = $cutoff;
  }

  if($self->{config}->{output_format} eq "vcf") {
    $output_vcf = 1;
  }

  if(defined($param_hash->{split_output})) {
    $self->{split_output} = 1;
  }

  return $self;
}

sub feature_types {
  return ["Transcript"];
}

sub get_header_info {
  my $self = shift;

  my %header;

  # Get the SpliceAI tool version from one of the files
  my $spliceai_version;
  my $spliceai_file = $self->{_files}[0];
  open(my $fh, "tabix -H $spliceai_file |") or die "Could not open '$spliceai_file': $!";
    while (my $line = <$fh>) {
      last if $line !~ /^#/;
      if ($line =~ /^##INFO/) {
          if ($line =~ /(SpliceAIv[\d.]+)/) {
            $spliceai_version = $1;
            $spliceai_version =~ s/SpliceAI//;
        }
      }
    }
  close($fh);

  if($output_vcf || $self->{split_output}) {
    $header{"SpliceAI_pred_SYMBOL"} = "SpliceAI ($spliceai_version) gene symbol";
    $header{"SpliceAI_pred_DS_AG"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta score for acceptor gain";
    $header{"SpliceAI_pred_DS_AL"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta score for acceptor loss";
    $header{"SpliceAI_pred_DS_DG"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta score for donor gain";
    $header{"SpliceAI_pred_DS_DL"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta score for donor loss";
    $header{"SpliceAI_pred_DP_AG"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta position for acceptor gain";
    $header{"SpliceAI_pred_DP_AL"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta position for acceptor loss";
    $header{"SpliceAI_pred_DP_DG"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta position for donor gain";
    $header{"SpliceAI_pred_DP_DL"} = "SpliceAI ($spliceai_version) predicted effect on splicing. Delta position for donor loss";
  }

  else {
    $header{"SpliceAI_pred"} = "SpliceAI ($spliceai_version) predicted effect on splicing. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL";
  }

  if($self->{cutoff}) {
    $header{"SpliceAI_cutoff"} = "Flag if delta score pass the cutoff (PASS) or if it does not (FAIL)";
  }

  return \%header;
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;
  my $chr = $vf->{chr};

  my $end = $vf->{end};
  my $start = $vf->{start};
  ($start, $end) = ($end, $start) if $start > $end;

  my @data = @{$self->get_data($chr, $start, $end)} if(defined $chr);

  return {} unless(@data);

  my $result_data = "";
  my $result_flag;

  # Store all SpliceAI results
  my %hash_aux;

  foreach my $data_value (@data) {

    my $ref_allele = $vf->ref_allele_string;
    my $allele_number = $tva->allele_number; # necessary for multi-allelic variants

    # Fetch the specific alt allele
    my @alt_alleles = $tva->base_variation_feature->alt_alleles;
    my $alt_allele = $alt_alleles[0][$allele_number - 1];

    # SpliceAI data is represented in VCF format
    # Convert VF ins/del to VCF to be able to compare alleles
    if($vf->allele_string =~ /-/) {
        my $convert_to_vcf = $vf->to_VCF_record;
        $start = ${$convert_to_vcf}[1];
        $ref_allele = ${$convert_to_vcf}[3];
        my $alt_allele_string = ${$convert_to_vcf}[4];
        my @aux_alt = split /,/, $alt_allele_string;
        $alt_allele = $aux_alt[$allele_number - 1];
    }

    if ($start == $data_value->{start} && $ref_allele eq $data_value->{ref} && $alt_allele eq $data_value->{alt}) {
      my %hash;

      if($output_vcf || $self->{config}->{output_format} eq "json" || $self->{config}->{rest} || $self->{split_output})  {
        my @data_values = split /\|/, $data_value->{result};
        my $prefix = "";
        $prefix = "SpliceAI_pred_" if($output_vcf || ($self->{split_output} && $self->{config}->{output_format} ne "json"));
        $hash{$prefix. "SYMBOL"} = $data_values[0];
        $hash{$prefix. "DS_AG"} = $data_values[1];
        $hash{$prefix. "DS_AL"} = $data_values[2];
        $hash{$prefix. "DS_DG"} = $data_values[3];
        $hash{$prefix. "DS_DL"} = $data_values[4];
        $hash{$prefix. "DP_AG"} = $data_values[5];
        $hash{$prefix. "DP_AL"} = $data_values[6];
        $hash{$prefix. "DP_DG"} = $data_values[7];
        $hash{$prefix. "DP_DL"} = $data_values[8];
      }

      else {
        $hash{"SpliceAI_pred"} = $data_value->{result};
      }

      # Add a flag if cutoff is used
      if($self->{cutoff}) {
        if($data_value->{info} >= $self->{cutoff}) {
          $result_flag = "PASS";
        }
        else {
          $result_flag = "FAIL";
        }
        $hash{"SpliceAI_cutoff"} = $result_flag;
      }

      $hash_aux{$data_value->{gene}} = \%hash;
    }
  }

  return {} unless(%hash_aux);

  my $result = {};

  # find the SpliceAI gene matching the variant gene symbol, if there is a match
  my $gene_symbol = $tva->transcript->{_gene_symbol} || $tva->transcript->{_gene_hgnc};
  if(($gene_symbol) && ($hash_aux{$gene_symbol})) {
      $result = ($self->{config}->{output_format} eq "json" || $self->{config}->{rest}) ?  {SpliceAI => $hash_aux{$gene_symbol}} : $hash_aux{$gene_symbol};
  }

  return $result;
}

# Parse data from SpliceAI file
sub parse_data {
  my ($self, $line) = @_;

  my ($chr, $start, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $line;

  $info =~ s/SpliceAI=//;
  my @info_splited = split (qr/\|/,$info, 3);
  my $allele = $info_splited[0];
  my $data = $info_splited[1] . "|" . $info_splited[2];
  my $gene = $info_splited[1];

  my $max_score;
  if($self->{cutoff}){
    my @scores = split (qr/\|/,$data);
    my @scores_list = @scores[1..4];
    $max_score = max(@scores_list);
  }

  return {
    chr    => $chr,
    start  => $start,
    ref    => $ref,
    alt    => $alt,
    info   => $max_score,
    result => $data,
    gene   => $gene,
  };
}

1;
