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

 Ensembl <https://www.ensembl.org/info/about/contact/index.html>
=cut

=head1 NAME

 G2P

=head1 SYNOPSIS

 mv G2P.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin G2P,file=/path/to/G2P.csv

=head1 DESCRIPTION

 A VEP plugin that uses G2P allelic requirements to assess variants in genes
 for potential phenotype involvement.

 The plugin has multiple configuration options, though minimally requires only
 the CSV file of G2P data.
 This Plugin is available for GRCh38 and GRCh37.

 For further information see:
 Thormann A, Halachev M, McLaren W, et al. Flexible and scalable diagnostic filtering of genomic variants using G2P with Ensembl VEP.
 Nature Communications. 2019 May;10(1):2373. doi:10.1038/s41467-019-10016-3. PMID: 31147538; PMCID: PMC6542828.


 Options are passed to the plugin as key=value pairs, (defaults in parentheses):

 file                  : Path to G2P data file. The file needs to be uncompressed.
                         - Download from https://www.ebi.ac.uk/gene2phenotype/downloads
                         - Download from PanelApp  

 variant_include_list  : A list of variants to include even if variants do not pass allele
                         frequency filtering. The include list needs to be a sorted, bgzipped and
                         tabixed VCF file.

 af_monoallelic        : maximum allele frequency for inclusion for monoallelic genes (0.0001)

 af_biallelic          : maximum allele frequency for inclusion for biallelic genes (0.005)
 confidence_levels     : Confidence levels include: definitive, strong, moderate, limited
                         Former confidence terms are still supported: confirmed, probable, possible, both RD and IF.
                         Separate multiple values with '&'.
                         https://www.ebi.ac.uk/gene2phenotype/terminology
                         Default levels are confirmed and probable.
 all_confidence_levels : Set to 1 to include all confidence levels
                         Setting the value to 1 will overwrite any confidence levels provided with the
                         confidence_levels option.
 af_from_vcf           : set value to 1 to include allele frequencies from VCF file. 
                         Specifiy the list of reference populations to include with '--af_from_vcf_keys'
 af_from_vcf_keys      : VCF collections used for annotating variant alleles with observed
                         allele frequencies. Allele frequencies are retrieved from VCF files. If
                         af_from_vcf is set to 1 but no VCF collections are specified with '--af_from_vcf_keys'
                         all available VCF collections are included. 
                         Available VCF collections: 'topmed', 'uk10k', 'gnomADe', 'gnomADe_r2.1.1', 'gnomADg', 'gnomADg_v3.1.2', 'gnomADev4.1', 'gnomADgv4.1'.
                         Separate multiple values with '&'.
                         VCF collections contain the following populations: 
                         * 'topmed' - TOPMed (available for GRCh37 and GRCh38).
                         * 'uk10k' - ALSPAC, TWINSUK (available for GRCh37 and GRCh38).
                         * 'gnomADe' & 'gnomADe_r2.1.1 & 'gnomADev4.1' - gnomADe:AFR, gnomADe:ALL, gnomADe:AMR, gnomADe:ASJ, gnomADe:EAS, gnomADe:FIN, gnomADe:NFE, gnomADe:OTH, gnomADe:SAS (for GRCh37 and GRCh38 respectively).
                         * 'gnomADg' & 'gnomADg_v3.1.2' & 'gnomADgv4.1' - gnomADg:AFR, gnomADg:ALL, gnomADg:AMR, gnomADg:ASJ, gnomADg:EAS, gnomADg:FIN, gnomADg:NFE, gnomADg:OTH (for GRCh37 and GRCh38 respectively).
                         Need to use 'af_from_vcf' parameter to use this option. 
  only_vcf_freq        : set to 1 to only use frequency from vcf files, can only be set if af_from_vcf is set.  
                         N/B - frequency information may be lost if this option is used 
 default_af            : default frequency of the input variant if no frequency data is
                         found (0). This determines whether such variants are included;
                         the value of 0 forces variants with no frequency data to be
                         included as this is considered equivalent to having a frequency
                         of 0. Set to 1 (or any value higher than 'af') to exclude them.
 types                 : SO consequence types to include. Separate multiple values with '&'
                         (splice_donor_variant, splice_acceptor_variant, stop_gained,
                         frameshift_variant, stop_lost, initiator_codon_variant,
                         inframe_insertion, inframe_deletion,missense_variant,
                         coding_sequence_variant, start_lost,transcript_ablation,
                         transcript_amplification, protein_altering_variant)
  
  log_dir              : write stats to log files in log_dir 

  txt_report           : write all G2P complete genes and attributes to txt file

  html_report          : write all G2P complete genes and attributes to html file

  filter_by_gene_symbol : set to 1 if filter by gene symbol.
                          Do not set if filtering by HGNC_id.
                          This option is set to 1 when using PanelApp files. 

  only_mane            : set to 1 to ignore transcripts that are not MANE
                         N/B - Information may be lost if this option is used.
               

 For more information - https://www.ebi.ac.uk/gene2phenotype/g2p_vep_plugin
 
 Example:
 --plugin G2P,file=G2P.csv,af_monoallelic=0.05,types=stop_gained&frameshift_variant
 --plugin G2P,file=G2P.csv,af_monoallelic=0.05,af_from_vcf=1
 --plugin G2P,file=G2P.csv,af_from_vcf=1,af_from_vcf_keys='topmed&gnomADe_r2.1.1'
 --plugin G2P,file=G2P.csv,af_from_vcf=1,af_from_vcf_keys='topmed&gnomADe_r2.1.1',confidence_levels='confirmed&probable&both RD and IF' 
 --plugin G2P,file=G2P.csv

=cut

package G2P;

use strict;
use warnings;
use Cwd;
use Scalar::Util qw(looks_like_number);
use FileHandle;
use Text::CSV;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::Variation::Utils::VEP qw(parse_line);
use Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);
use List::Util qw(any);

our $CAN_USE_HTS_PM;

BEGIN {
  if (eval { require Bio::DB::HTS::Tabix; 1 }) {
    $CAN_USE_HTS_PM = 1;
  }
}


my %DEFAULTS = (

  # vars must have a frequency <= to this to pass
  af_monoallelic => 0.0001,
  af_biallelic => 0.005, 

  af_keys => [qw(AA AFR AMR EA EAS EUR SAS gnomAD gnomAD_AFR gnomAD_AMR gnomAD_ASJ gnomAD_EAS gnomAD_FIN gnomAD_NFE gnomAD_OTH gnomAD_SAS
                 gnomADg gnomADg_AFR gnomADg_AMR gnomADg_ASJ gnomADg_EAS gnomADg_FIN gnomADg_NFE gnomADg_OTH gnomADg_SAS gnomADg_EAS gnomADg_REMAINING gnomADg_MID
                 gnomADe gnomADe_AFR gnomADe_AMR gnomADe_ASJ gnomADe_EAS gnomADe_FIN gnomADe_NFE gnomADe_OTH gnomADe_SAS gnomADe_EAS gnomADe_REMAINING gnomADe_MID)],

  af_from_vcf_keys => [qw(uk10k topmed gnomADe gnomADe_r2.1.1 gnomADg gnomADg_v3.1.2 gnomADe_v4.1 gnomADg_v4.1)],

  # if no MAF data is found, default to 0
  # this means absence of MAF data is considered equivalent to MAF=0
  # set to 1 to do the "opposite", i.e. exclude variants with no MAF data
  default_af => 0,
  # adding new confidence levels based on the new terminology 
  confidence_levels => [qw(confirmed probable definitive strong moderate)],

  # only include variants with these consequence types
  # currently not ontology-resolved, exact term matches only
  types => {map {$_ => 1} qw(splice_donor_variant splice_acceptor_variant stop_gained frameshift_variant stop_lost initiator_codon_variant inframe_insertion inframe_deletion missense_variant coding_sequence_variant start_lost transcript_ablation transcript_amplification protein_altering_variant)},

);

my $af_key_2_population_name = {
  minor_allele_freq => 'global allele frequency (AF) from 1000 Genomes Phase 3 data',
  AFR => '1000GENOMES:phase_3:AFR',
  AMR => '1000GENOMES:phase_3:AMR',
  EAS => '1000GENOMES:phase_3:EAS',
  EUR => '1000GENOMES:phase_3:EUR',
  SAS => '1000GENOMES:phase_3:SAS',
  AA => 'Exome Sequencing Project 6500:African_American',
  EA => 'Exome Sequencing Project 6500:European_American',
  gnomAD => 'Genome Aggregation Database:Total',
  gnomAD_AFR => 'Genome Aggregation Database exomes v2.1:African/African American',
  gnomAD_AMR => 'Genome Aggregation Database exomes v2.1:Latino/Admixed American',
  gnomAD_ASJ => 'Genome Aggregation Database exomes v2.1:Ashkenazi Jewish',
  gnomAD_EAS => 'Genome Aggregation Database exomes v2.1:East Asian',
  gnomAD_FIN => 'Genome Aggregation Database exomes v2.1:Finnish',
  gnomAD_NFE => 'Genome Aggregation Database exomes v2.1:Non-Finnish European',
  gnomAD_OTH => 'Genome Aggregation Database exomes v2.1:Other (population not assigned)',
  gnomAD_SAS => 'Genome Aggregation Database exomes v2.1:South Asian',
  gnomADe_MID => "Genome Aggregation Database exomes Mid Eastern",
  gnomADg_MID => "Genome Aggregation Database genomes Mid Eastern",
  gnomADe_REMAINING => "Genome Aggregation Database exomes Remaining",
  gnomADg_REMAINING => "Genome Aggregation Database genomes Remaining",
};

my $allelic_requirements = {
  'biallelic' => { af => 0.005, rules => {HET => 2, HOM => 1} },
  'biallelic_autosomal' => { af => 0.005, rules => {HET => 2, HOM => 1} },
  'biallelic_PAR' => { af => 0.005, rules => {HET => 2, HOM => 1} },
  'monoallelic' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'monoallelic_autosomal' =>  => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'hemizygous' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'hemizygous_biallelic' => { af => 0.0001, rules => {HET => 2, HOM => 1} },
  'monoallelic_X_hem'  => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'monoallelic_Y_hem'  => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'x-linked dominant' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'monoallelic_X_het' =>  { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'x-linked over-dominance' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'mitochondrial' =>  { af => 0.0001, rules => {HET => 1, HOM => 1} },
  'monoallelic_PAR' => { af => 0.0001, rules => {HET => 1, HOM => 1} },
};

my $supported_confidence_levels = {
  'confirmed' => 1,
  'definitive' => 1,
  'probable' => 1,
  'strong' => 1,
  'moderate' => 1, 
  'possible' => 1,
  'limited' => 1,
  'both RD and IF' => 1,
};

my @allelic_requirement_terms = keys %$allelic_requirements;

# keys containing the assembly and the key to do a quick key and assembly lookup
my $afvcf_keys = {
    "uk10k_GRCh37" => 1, 
    "uk10k_GRCh38" => 1, 
    "topmed_GRCh37" => 1, 
    "topmed_GRCh38" => 1,
    "gnomADe_GRCh37" => 1, 
    "gnomADe_r2.1.1_GRCh38" => 1,
    "gnomADg_GRCh37" => 1, 
    "gnomADg_v3.1.2_GRCh38" => 1,
    "gnomADg_v4.1_GRCh37" => 1,
    "gnomADe_v4.1_GRCh37" => 1,
    "gnomADg_v4.1_GRCh38" => 1,
    "gnomADe_v4.1_GRCh38" => 1,
};

my $a_keys = {
  "uk10k" => 1,
  "topmed" => 1,
  "gnomADe" => 1,
  "gnomADe_r2.1.1" => 1,
  "gnomADg" => 1,
  "gnomADg_v3.1.2" => 1,
  "gnomADg_v4.1" => 1,
  "gnomADe_v4.1" => 1
};

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
# suppress warnings that the FeatureAdpators spit if using no_slice_cache
  Bio::EnsEMBL::Utils::Exception::verbose(1999);

  my $params = $self->params_to_hash();
  my $file = '';


  # user only supplied file as first param?
  if (!keys %$params) {
    $file = $self->params->[0];
    $params->{file} = $file;
  }
  else {
    $file = $params->{file};
    open IN, "<",  $file;
    # supporting panelapp by always filtering by gene symbol if PanelApp file 
    while (<IN>){
      $params->{filter_by_gene_symbol} = 1 if (/Model_Of_Inheritance/);
      last;
    }
    close $file;

    # process types
    if ($params->{types}) {
      $params->{types} = {map {$_ => 1} split(/[\;\&\|]/, $params->{types})};
    }

    # check af
    foreach my $af (qw/af_monoallelic af_biallelic/) {
      if($params->{$af}) {
        die("ERROR: Invalid value for af: ".$params->{$af} . "\n") unless
          looks_like_number($params->{$af}) && ($params->{$af} >= 0 && $params->{$af} <= 1)
      }
      my $ar = $af;
      $ar =~ s/af_//;
      $allelic_requirements->{$ar}->{af} = $params->{$af} if (defined $params->{$af});
    }

    $params->{af_keys} = \@{$DEFAULTS{af_keys}};
  }

  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
  $year += 1900;
  $mon++;
  my $stamp = join('_', ($year, $mon, $mday, $hour, $min, $sec));
  my $cwd_dir = getcwd;
  my $new_log_dir = "$cwd_dir/g2p_log_dir\_$stamp";
  my $log_dir = $params->{log_dir} || $new_log_dir;
  if (!-d $log_dir) {
    my $return = mkdir $log_dir, 0755;
    die("ERROR: Couldn't create log_dir $log_dir $!\n") if (!$return);
    $params->{log_dir} = $log_dir;
  } 
  else{
    opendir my $dh, $log_dir or die("ERROR: There was a problem opening the log_dir: $!\n");
    my @check = grep {$_ ne '.' and $_ ne '..'} readdir $dh;
    die("ERROR: The log directory ($log_dir) is not empty. You need to empty directory before using the plugin \n") if (scalar @check != 0);
    closedir $dh;
    $params->{log_dir} = $log_dir;
  }

  foreach my $report_type (qw/txt_report html_report/) {
    if (!$params->{$report_type}) {
      my $file_type = ($report_type eq 'txt_report') ? 'txt' : 'html';
      $params->{$report_type} = $cwd_dir . "/$report_type\_$stamp.$file_type";
    } 
  }

  if ($params->{all_confidence_levels}) {
    if ($params->{confidence_levels}) {
      warn("Option all_confidence_levels set to 1 overwrites confidence levels provided with confidence_levels option.");
    }
    $params->{confidence_levels} = ['possible', @{$DEFAULTS{confidence_levels}}];
  }
  elsif ($params->{confidence_levels}) {
    my @confidence_levels = ();
    foreach my $confidence_level (split(/[\;\&\|]/, $params->{confidence_levels})) {
      if (!$supported_confidence_levels->{$confidence_level}) {
        die "$confidence_level is not a supported value for supported confidence levels. Supported values are: ", join(', ', keys %$supported_confidence_levels);
      } else {
        push @confidence_levels, $confidence_level;
        push @confidence_levels, 'both DD and IF' if ($confidence_level eq 'both RD and IF'); # legacy support for using both DD and IF

      }
    }
    if (scalar @confidence_levels > 0) {
      $params->{confidence_levels} = \@confidence_levels;
    }
  }
  
  die "only_vcf_freq option can only be used with af_from_vcf"  if ($params->{only_vcf_freq} && !$params->{af_from_vcf});
  if ($params->{af_from_vcf}) {
    if ($CAN_USE_HTS_PM) {
      my @vcf_collection_ids = ();
      # adding a die if assembly is not used and af_from_vcf keys option is used
      my $assembly =  $self->{config}->{assembly};
      die "Assembly needs to be defined to use af_from_vcf option" if (!defined ($assembly));
      if ($params->{af_from_vcf_keys}) {
        foreach my $key (split(/[\;\&\|]/, $params->{af_from_vcf_keys})) {
          my $key_assembly = $key."_".$assembly;
          if (!$afvcf_keys->{$key_assembly}){
            # to die if key is not supported, checking with the key and the assembly 
            die "$key is not a supported key. Supported keys and assembly are: ", join(',', keys %$a_keys), ".\n
                 gnomADe, gnomADg, gnomADg_v4.1 and gnomADe_v4.1 is supported for assembly GRCh37 \n
                 gnomADe_r2.1.1, gnomADg_v3.1.2 and gnomADg_v4.1 and gnomADe_v4.1 is supported for assembly GRCh38 \n" ;
          }
          else {
            push @vcf_collection_ids, $key;
            push @vcf_collection_ids, $key_assembly;
          }
        }
      } else {
        foreach my $key (@{$DEFAULTS{af_from_vcf_keys}}) {
          push @vcf_collection_ids, $key;
          push @vcf_collection_ids, "$key\_$assembly";
        }
      }

      my $species =  $self->{config}->{species};
      my $reg = $self->{config}->{reg};
      my $vca;
      if (defined $self->{config}->{offline}) {
        $vca = Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor->new();
      } else {
        my $vdba = $reg->get_DBAdaptor($species, 'variation');
        $vdba->dbc->reconnect_when_lost(1);
        $vca = $vdba->get_VCFCollectionAdaptor;
        $vca->db->use_vcf(2);
      }
      my $vcf_collections = $vca->fetch_all;
      my @collections = ();
      foreach my $vcf_collection (@$vcf_collections) {
        $vcf_collection->use_db(0) if (defined $self->{config}->{offline});
        my $vcf_collection_id = $vcf_collection->id;
        if ($vcf_collection->assembly eq $assembly && grep {$_ =~ /$vcf_collection_id/i} @vcf_collection_ids) {
          delete $vcf_collection->adaptor->{collections};
          delete $vcf_collection->adaptor->{config};
          my $description = $vcf_collection->description || $vcf_collection_id;
          foreach my $population (@{$vcf_collection->get_all_Populations})  {
            my $population_name = $population->name;
            my $population_description = $population->description;    
            $af_key_2_population_name->{$population_name} = "$description $population_name $population_description";
          }
          push @collections, $vcf_collection;
        }
      }
      warn "Couldn't find VCF collection ids for assembly " . $assembly if (!@collections);
      $self->{config}->{vcf_collections} = \@collections;
      $self->{config}->{use_vcf} = 1;
      $self->{config}->{use_only_vcf} = 1 if $params->{only_vcf_freq};
    } else {
      warn "Cannot get data from VCF without Bio::DB::HTS::Tabix";
    } 
  }

  if ($params->{variant_include_list}) {
    if (! -f $params->{variant_include_list}) {
      die "Variant include list (" . $params->{variant_include_list} . ") does not exist.";
    }
    $self->{_files} = [$params->{variant_include_list}];
  }
  
  if (defined($params->{filter_by_gene_symbol}) and $params->{filter_by_gene_symbol} != 1) {
    $params->{filter_by_gene_symbol} = undef;
    die "The option --filter_by_gene_symbol needs to be set to 1 \n";
  }
  
  if (defined($params->{only_mane}) and $params->{only_mane} != 1) {
     $params->{only_mane} = undef;
    die "The option only_mane needs to be set to 1 \n";
  }
  
  if (defined($params->{only_mane}) and $self->{config}->{assembly} ne "GRCh38") {
    die "The option only_mane only works with GRCh38 assembly \n";
  }


  # copy in default params
  $params->{$_} //= $DEFAULTS{$_} for keys %DEFAULTS;
  $self->{user_params} = $params;

  $self->{config}->{frequency_threshold} = _get_highest_frequency_threshold();
  # read data from file
  $self->{gene_data} = $self->read_gene_data_from_file($file);
  $self->synonym_mappings();
  $self->hgnc_mappings();

  # force some config params
  $self->{config}->{individual} //= ['all'];
  $self->{config}->{symbol} = 1;
  $self->{config}->{check_existing} = 1;
  $self->{config}->{failed} = 1;
  $self->{config}->{af} = 1;
  $self->{config}->{af_1kg} = 1;
  $self->{config}->{af_gnomad} = 1;
  $self->{config}->{sift} = 'b';
  $self->{config}->{polyphen} = 'b';

  # tell VEP we have a cache so stuff gets shared/merged between forks
  $self->{has_cache} = 1;
  $self->{cache}->{g2p_in_vcf} = {};

  return $self;
}

=head2 _get_highest_frequency_threshold

  Description: Retrieve the highest allele frequency threshold across all defined allelic requirements.
               This will speed up the filtering by allele frequency. As soon as we found a frequency
               higher than this highest frequency the filtering fails and we don't need to consider
               the variant further.
  Returntype : Float $highest_frequency
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _get_highest_frequency_threshold {
  my $highest_frequency = 0.0;
  foreach my $ar (keys %$allelic_requirements) {
    if ($allelic_requirements->{$ar}->{af} > $highest_frequency) {
      $highest_frequency = $allelic_requirements->{$ar}->{af};
    }
  }
  return $highest_frequency;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  return {
    G2P_flag => 'Flags zygosity of valid variants for a G2P gene',
    G2P_complete => 'Indicates this variant completes the allelic requirements for a G2P gene',
    G2P_gene_req => 'MONO or BI depending on the context in which this gene has been explored',
  };
}

=head2 run

  Arg [1]    : TranscriptVariationAllele $tva
  Arg [2]    : Hashref $line
  Description: Filter input transcript variation allele on:
                - G2P gene overlap
                - variant consequence
                - variant include list
                - allele frequency
               Based on the filtering results check if the allelic requirement is fulfilled and write results to hash.
               Dump annnotations to a log file for generating TXT and HTML output files after VEP has finished.
  Returntype : Hashref $results
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub run {
  my ($self, $tva, $line) = @_;

  # only interested if we know the zygosity
  my $zyg = defined($line->{Extra}) ? $line->{Extra}->{ZYG} : $line->{ZYG};
  return {} unless $zyg;

  return {} if ($self->{user_params}->{only_mane} && !$tva->transcript->is_mane);
  # filter by G2P gene overlap
  return {} if (!$self->gene_overlap_filtering($tva));

  $self->set_variant_include_list_flag($tva);

  # filter by variant consequence
  return {} if (!$self->consequence_filtering($tva));
  
  # filter by allele frequency
  return {} if (!$self->frequency_filtering($tva));

  # dump annotations for txt and html report files
  $self->dump_vf_annotations($tva);      
  $self->dump_individual_annotations($tva, $zyg);

  # check if transcript contains enough variants to fulfill the allelic requirement of the gene
  my $G2P_complete = $self->is_g2p_complete($tva, $zyg);
  my $G2P_flag = $self->is_valid_g2p_variant($tva, $zyg);
  my $results = {};
  $results->{G2P_complete} = $G2P_complete if ($G2P_complete); 
  $results->{G2P_flag} = $G2P_flag if ($G2P_flag);
  return $results;
}

=head2 set_variant_include_list_flag

  Arg [1]    : TranscriptVariationAllele $tva
  Description: Check if variant is part of the variant include list.
               If the variant is on the variant include list then report
               it in the result set regardless if the variant passed the
               filtering by variant consequence and allele frequency.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub set_variant_include_list_flag {
  my $self = shift;
  my $tva = shift;
  return if (!$self->{user_params}->{variant_include_list});
  my $vf = $tva->variation_feature;

  my $allele = $tva->base_variation_feature->alt_alleles;

  foreach (@{$self->get_data($vf->{chr}, $vf->{start} - 1, $vf->{end})}) {
    my @vcf_alleles = split /\//, $_->allele_string;
    my $ref_allele  = shift @vcf_alleles;
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $allele,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
        ref  => $ref_allele,
        alts => \@vcf_alleles,
        pos  => $_->{start},
      }
    );
    if (scalar @$matches) {
      my $vf_cache_name = $self->get_cache_name($vf);
      $self->{g2p_vf_cache}->{$vf_cache_name}->{is_on_variant_include_list} = 1;
      last;
    }
  }
}

=head2 is_valid_g2p_variant

  Arg [1]    : TranscriptVariationAllele $tva
  Arg [2]    : String $zygosity
  Example    : $valid_g2p_variant = $self->is_valid_g2p_variant($tva, 'HOM')
  Description: Take all allelic requirements of the gene that overlap this variant
               and check if the variant passes the frequency threshold filter where the threshold is defined
               by the allelic requirement of the gene.
               Concatenate results for several allelic requirements by ','.
  Returntype : String for example "monoallelic=HOM"
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub is_valid_g2p_variant {
  my $self = shift;
  my $tva = shift;
  my $zyg = shift;
  my $transcript = $tva->transcript;
  my $gene_stable_id = get_gene_stable_id($transcript);
  my @allelic_requirements = keys %{$self->{ar}->{$gene_stable_id}};
  my @results = ();
  foreach my $ar (@allelic_requirements) {
    my $ar_rules = $allelic_requirements->{$ar};
    my $af_threshold = $ar_rules->{af};
    my $variants = $self->variants_filtered_by_frequency_threshold($af_threshold, [$self->{vf_cache_name}]);
    if (scalar @$variants > 0) {
      push @results, "$ar=$zyg";
    }
  }
  return join(',', @results);
}

=head2 is_g2p_complete

  Arg [1]    : TranscriptVariationAllele $tva
  Arg [2]    : String $zygosity           
  Description: A G2P gene is considered complete if its allelic requirement is fulfilled.
               Create a summary string which connects the fulfilled allelic_requirement and all variants which pass filtering. 
  Returntype : String $g2p_complete, for example "monoallelic=HET:7_941481_C/T&HET:7_929274_A/G,HOM:7_931481_C/T|biallelic=HOM:7_931481_C/T"
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub is_g2p_complete {
  my $self = shift;
  my $tva = shift;
  my $zyg = shift;
  my $vf = $tva->base_variation_feature;
  my $individual = $vf->{individual};
  my $transcript = $tva->transcript;
  my $gene_stable_id = get_gene_stable_id($transcript);
  my $transcript_stable_id = $transcript->stable_id; 
  $self->{per_individual}->{$individual}->{$transcript_stable_id}->{$zyg}->{$self->{vf_cache_name}} = 1;
  my @allelic_requirements = keys %{$self->{ar}->{$gene_stable_id}};
  my @g2p_complete = ();
  foreach my $ar (@allelic_requirements) {
    my $zyg2var = $self->{per_individual}->{$individual}->{$transcript_stable_id};
    my $filtered_zyg2var = $self->zyg2var_filtered_by_allelic_requirement_rule($ar, $zyg2var);
    if (defined $filtered_zyg2var) {
      my @filtered_variants = ();
      foreach my $zyg (keys %$filtered_zyg2var) {
        push @filtered_variants, join('&', map {"$zyg:$_"} @{$filtered_zyg2var->{$zyg}});
      }
      push @g2p_complete, "$ar=" . join(',', @filtered_variants);
    }
  }
  return join('\|', @g2p_complete);
} 

=head2 zyg2var_filtered_by_allelic_requirement_rule

  Arg [1]    : String $ar allelic requirement
  Arg [2]    : Hashref $zyg2var
  Example    : $zyg2var_filtered = $self->zyg2var_filtered_by_allelic_requirement_rule('biallelic', {
                'HET' => {
                  '4_32941481_C/T' => 1,
                  '4_32929274_A/G' => 1
                }
               });
  Description: Check if the variants fulfil the given allelic requirement. Variants are grouped by their
               zygosity. The subroutine checks for each variant if the internally stored frequency for a variant
               is lower than the allele frequency threshold defined by the allelic requirement. Then consider
               the variants which pass frequency filtering and check if the number is sufficient to fulfil the
               allelic requirement.  
  Returntype : Hashref $zyg2var_filtered: with zygosity as key and filtered variants in arrayref as value
               For example: {
                'HET' => [
                  '4_32941481_C/T',
                  '4_32929274_A/G'
                ]
               };
               Undef, if none of the variants pass the filtering
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub zyg2var_filtered_by_allelic_requirement_rule {
  my $self = shift;
  my $ar = shift;
  my $zyg2variants = shift;
  my $ar_rules = $allelic_requirements->{$ar};
  # allele frequency threshold for given allelic requirement
  my $af_threshold = $ar_rules->{af};
  # number of variants with a particular zygosity that need to pass filtering
  # to fulfil allelic requirement
  my $zyg2counts = $ar_rules->{rules};
  my $results = {};
  foreach my $zyg (keys %$zyg2counts) {
    my $count = $zyg2counts->{$zyg};
    my @all_variants = keys %{$zyg2variants->{$zyg}};
    my $variants = $self->variants_filtered_by_frequency_threshold($af_threshold, \@all_variants);
    if (scalar @$variants >= $count) {
      $results->{$zyg} = $variants;
    }
  }
  if (scalar keys %$results > 0) {
    return $results;
  } else {
    return undef;
  }
}

=head2 variants_filtered_by_frequency_threshold

  Arg [1]    : Float $af_threshold allele frequency threshold
  Arg [2]    : Arrayref $variants
  Example    : $variants_filtered = $self->variants_filtered_by_frequency_threshold(0.001, [
                  '4_32941481_C/T' => 1,
                  '4_32929274_A/G' => 1
                ]
               );
  Description: Check for each variant if the highest frequency that has been observed in any
               population (stored internally under the highest_frequencies hash key for a variant) is lower
               than the given allele frequency threshold. If yes, add the variant to the result set.
               Also add the variant to the result set if the variant hasn't any observed allele frequencies
               or if the variant is in the variant include list. 
  Returntype : Arrayref $variants_filtered
               For example: [
                  '4_32941481_C/T' => 1,
                  '4_32929274_A/G' => 1
               ];
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub variants_filtered_by_frequency_threshold {
  my $self = shift;
  my $af_threshold = shift;
  my $variants = shift;
  my @pass_variants = ();
  foreach my $variant (@$variants) {
    # get the highest frequency that has been observed for the variant
    my $highest_frequency = $self->highest_frequency($variant);
    if (! defined $highest_frequency || (defined $af_threshold && $highest_frequency <= $af_threshold) ||
         $self->{g2p_vf_cache}->{$variant}->{is_on_variant_include_list}
    ) {
      push @pass_variants, $variant;
    }
  }
  return \@pass_variants;
}

=head2 gene_overlap_filtering

  Arg [1]    : TranscriptVariationAllele $tva
  Example    :
  Description: returns 1 or 0 depending on if the gene is part of the input panel. If the gene is part of the panel we store
               the allelic requirement in the internal cash under the name 'ar'. We also write G2P_gene_data and G2P_in_vcf
               information to the log file. We call _dump_transcript_annotations and write transcript information to the log file.
  Returntype : Boolean
  Exceptions : None
  Caller     : run
  Status     : Stable

=cut
sub gene_overlap_filtering {
  my $self = shift;
  my $tva = shift;
  my $transcript = $tva->transcript;
  my $gene_stable_id = get_gene_stable_id($transcript);
  my $pass_gene_overlap_filter = $self->{g2p_gene_cache}->{$gene_stable_id};
  my @gene_xrefs = ();
  if (! defined $pass_gene_overlap_filter && ! defined $self->{user_params}->{filter_by_gene_symbol}) {
    my @gene_hgnc = split /:/, $transcript->{_gene_hgnc_id} if (defined $transcript->{_gene_hgnc_id});
    my $hgnc_id = $gene_hgnc[1] if (@gene_hgnc);
    if (defined $hgnc_id) {
      foreach my $id (keys %{ $self->{hgnc_mapping} }) {
        if ( $hgnc_id == $id) {
          my $gene_symbol = $self->{hgnc_mapping}{$id};
          my $gene_data = $self->gene_data($gene_symbol) if defined ($gene_symbol);
          if (defined $gene_data) {
            if (defined $gene_data->{'allelic requirement'} && scalar @{$gene_data->{'allelic requirement'}}) {
              foreach my $ar (@{$gene_data->{'allelic requirement'}}) {
                $self->{ar}->{$gene_stable_id}->{$ar} = 1;
              }
              $self->write_report('G2P_gene_data', $gene_stable_id, $gene_data, $gene_data->{'gene_xrefs'}, $gene_data->{'HGNC'}, $gene_data->{'confidence_category'}, $gene_data->{'confidence_value'} );
            }
            $self->write_report('G2P_in_vcf', $gene_stable_id);
            $pass_gene_overlap_filter = 1;
            last;
          }
        } 
      }
    }
    $self->{g2p_gene_cache}->{$gene_stable_id} = $pass_gene_overlap_filter;
  }
  if (! defined $pass_gene_overlap_filter && defined $self->{user_params}->{filter_by_gene_symbol} == 1) {
    my $gene_symbol = $transcript->{_gene_symbol} || $transcript->{_gene_hgnc};
    $pass_gene_overlap_filter = 0;
    foreach my $gene_id ($gene_symbol, $gene_stable_id) {
      my $gene_data = $self->gene_data($gene_id) if (defined $gene_id);
      if (defined $gene_data) {
        if (defined $gene_data->{'allelic requirement'} && scalar @{$gene_data->{'allelic requirement'}}) {
          foreach my $ar (@{$gene_data->{'allelic requirement'}}) {
            $self->{ar}->{$gene_stable_id}->{$ar} = 1;
          } 
          $self->write_report('G2P_gene_data', $gene_stable_id, $gene_data, $gene_data->{'gene_xrefs'}, undef,  $gene_data->{'confidence_category'}, $gene_data->{'confidence_value'});
        } 
        $self->write_report('G2P_in_vcf', $gene_stable_id);
        $pass_gene_overlap_filter = 1;
        last;
      } 
    }
    $self->{g2p_gene_cache}->{$gene_stable_id} = $pass_gene_overlap_filter;
  }
  $self->_dump_transcript_annotations($transcript) if ($pass_gene_overlap_filter);
  return $self->{g2p_gene_cache}->{$gene_stable_id};
}

=head2 consequence_filtering

  Arg [1]    : TranscriptVariationAllele $tva
  Description: returns 1 or 0 depending on if the variant passes consequence filtering.
               If the variant is on the variant include list we return 1, regardless
               if the variant consequence is defined in types.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub consequence_filtering {
  my $self = shift;
  my $tva = shift;
  my $vf = $tva->base_variation_feature;
  my $vf_cache_name = $self->get_cache_name($vf);
  return ((grep { $self->{user_params}->{types}->{$_->SO_term} } @{$tva->get_all_OverlapConsequences}) ||
          $self->{g2p_vf_cache}->{$vf_cache_name}->{is_on_variant_include_list});
}

=head2 _dump_transcript_annotations

  Arg [1]    : Transcript $transcript
  Description: Write G2P_transcript_data to log file.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _dump_transcript_annotations {
  my $self = shift;
  my $transcript = shift;
  my $transcript_stable_id = $transcript->stable_id;
  if (!defined $self->{g2p_transcript_cache}->{$transcript_stable_id}) {
    my $gene_stable_id = get_gene_stable_id($transcript);
    if ($transcript->is_canonical) {
      $self->write_report('G2P_transcript_data', "$gene_stable_id\t$transcript_stable_id\tis_canonical");
    }
    $self->{g2p_transcript_cache}->{$transcript_stable_id} = 1;
  }
}

=head2 get_cache_name

  Arg [1]    : VariationFeature $vf
  Description: Create a variant identifier for internal use which is stored in the internal cache
  Returntype : String $vf_cache_name for example 19_8412862_G/A 
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub get_cache_name {
  my $self = shift;
  my $vf = shift;
  my $cache_name = ($vf->{original_chr} || $vf->{chr}) . '_' . $vf->{start} . '_' . ($vf->{allele_string} || $vf->{class_SO_term});
  return $cache_name;
}

=head2 get_gene_stable_id

  Arg [1]    : Transcript $transcript
  Description: Retrives the Ensembl or RefSeq gene stable id from the transcript object
  Returntype : String $gene_stable_id
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub get_gene_stable_id {
  my $transcript = shift;
  return $transcript->{_gene}->stable_id;
  
}

=head2 frequency_filtering

  Arg [1]    : TranscriptVariationAllele $tva
  Description: returns 1 or 0 depending on if the variant passes frequency filtering.
               We consider allele frequencies from the VEP cache and VCF files for filtering.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub frequency_filtering {
  my $self = shift;
  my $tva = shift;

  my $vf = $tva->base_variation_feature;
  # Set up caching to avoid looking up frequencies for each overlapping transcript
  my $vf_cache_name = $self->get_cache_name($vf);
  $self->{vf_cache_name} = $vf_cache_name;
  $self->{g2p_vf_cache} = {} if (!defined $self->{g2p_vf_cache}->{$vf_cache_name});
  # Retrieve cached result
  my $pass_frequency_filter = $self->{g2p_vf_cache}->{$vf_cache_name}->{pass_frequency_filter};
  return $pass_frequency_filter if (defined $pass_frequency_filter);
  # Check frequencies from cache files first
  $pass_frequency_filter = $self->_vep_cache_frequency_filtering($tva) if (!defined($self->{config}->{use_only_vcf}));
  # Check frequencies from VCF files if user is providing use_vcf flag
  if ($self->{config}->{use_vcf}) {
    $pass_frequency_filter = $self->_vcf_frequency_filtering($tva);
  } 

  $self->{g2p_vf_cache}->{$vf_cache_name}->{pass_frequency_filter} = $pass_frequency_filter;
  return $self->{g2p_vf_cache}->{$vf_cache_name}->{pass_frequency_filter};
}

=head2 _vep_cache_frequency_filtering

  Arg [1]    : TranscriptVariationAllele $tva
  Description: Returns 1 or 0 depending on if the variant passes frequency filtering where allele frequencies come from the VEP cache.
               Return 1 if no observed allele frequencies exist for the given variant which means that the filtering passes.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _vep_cache_frequency_filtering {
  my $self = shift;
  my $tva = shift;

  my $allele = $tva->base_variation_feature->alt_alleles;
  my $vf     = $tva->base_variation_feature;
  my $frequency_threshold = $self->{config}->{frequency_threshold}; 
  my $existing = $vf->{existing}; # Get existing variants from cache file which are stored on VF level
  my @keys = @{$self->{user_params}->{af_keys}}; # Consider user defined list of af keys
  my $dumped_annotations = 0;  # Indicates if existing annotations have already been dumped for txt and html report files
  my $vf_cache_name =  $self->{vf_cache_name};
  foreach my $existing_var (@$existing) {
    my @frequencies = grep defined, @{$existing_var}{@keys};
    if ($existing_var->{matched_alleles}) { # Get matched alleles from input variant and existing variant, in case input variant was normalized to match variant from cache file
      $allele = $existing_var->{matched_alleles}[0]->{b_allele};
    }

    next if (!@frequencies);
    if ($self->_exceeds_frequency_threshold(\@frequencies, $allele, $frequency_threshold) && !$self->{g2p_vf_cache}->{$vf_cache_name}->{is_on_variant_include_list}) { 
      return 0; # Return 0 (failed filtering) if frequencies exceed threshold and variant is not on variant_include_list
    } else {
      # Dump annotations for txt and html report files
      $self->_dump_existing_vf_frequencies($existing_var, $allele);
      $self->_dump_existing_vf_annotations($existing_var);
      $dumped_annotations = 1;
    }
  }
  # If we get to this point it means that there were no frequencies for the input variant in the cache files
  # and we pass the filtering step.
  # We need to dump 'empty' annotations for such variants to indicate that there are no available frequencies
  $self->_dump_existing_vf_annotations() if (!$dumped_annotations);
  return 1;
}

=head2 _dump_existing_vf_frequencies

  Arg [1]    : Hashref $existing_var
  Description: The $existing_var contains everything that is stored for the variant in the VEP cache file.
               Extract allele frequencies from $existing_var and write it to the log file as G2P_frequencies. 
               Store the highest observed frequency for the variant which is used later for filtering.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _dump_existing_vf_frequencies {
  my $self = shift;
  my $existing_var = shift;
  my $allele = shift;
  my @keys = @{$self->{user_params}->{af_keys}};
  my @frequencies = ();
  my $higest_frequency = 0;
  foreach my $population_name (@keys) {
    my $af = $existing_var->{$population_name};
    next if (!defined $af);
    foreach my $pair (split(',', $af)) {
      my ($a, $f) = split(':', $pair);
      if(($a || '') eq $allele && defined($f)) {
        push @frequencies, "$population_name=$f";
        $higest_frequency = $f if ($f > $higest_frequency);
      }
    }
  }
  $self->highest_frequency($self->{vf_cache_name}, $higest_frequency);
  $self->write_report('G2P_frequencies', $self->{vf_cache_name}, \@frequencies);
}

=head2 _dump_exisiting_vf_annotations

  Arg [1]    : Hashref $existing_var
  Description: The $existing_var contains everything that is stored for the variant in the VEP cache file.
               Extract variant annotation from the hashref and write the annotations to the log file as G2P_existing_vf_annotations. 
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _dump_existing_vf_annotations {
  my $self = shift;
  my $existing_var = shift;

  my $data = {
    'clin_sig' => 'NA',
    'failed' => 'NA',
    'existing_name' => 'NA',
    'novel' => 'yes',
  };
  if ($existing_var) { 
    $data = {
      'clin_sig' => $existing_var->{clin_sig} || 'NA',
      'failed' => ($existing_var->{failed}) ? 'yes' : 'no',
      'existing_name' => $existing_var->{variation_name} || 'NA',
      'novel' => 'no',
    };
  }
  $self->write_report('G2P_existing_vf_annotations', $self->{vf_cache_name}, $data);
}

=head2 _exceeds_frequency_threshold

  Arg [1]    : Arrayref $vep_cache_frequencies
  Arg [2]    : String $allele
  Arg [3]    : Float $threshold
  Example    : $exceeds_threshold = $self->_exceeds_frequency_threshold(['A:0.001', 'A:0.0001'], 'A', 0.01);
  Description: Returns 1 if any of the vep cache frequencies exceeds the given frequency threshold for the given allele. 
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _exceeds_frequency_threshold {
  my $self = shift;
  my $vep_cache_frequencies = shift;
  my $allele = shift;
  my $threshold = shift;
  foreach my $vep_cache_frequency (@$vep_cache_frequencies) {
    foreach my $pair (split(',', $vep_cache_frequency)) {
      my ($a, $f) = split(':', $pair);
      if(($a || '') eq $allele && defined($f)) {
        return 1 if ($f > $threshold);
      }
    }
  }
  return 0;
}

=head2 _vcf_frequency_filtering

  Arg [1]    : TranscriptVariationAllele $tva
  Description: Returns 1 or 0 depending on if the variant passes frequency filtering where allele frequencies come from VCF files.
               Return 1 if no observed allele frequencies exist for the given variant which means that the filtering passes.
  Returntype : Boolean
  Exceptions : None
  Caller     : General
  Status     : Stable


=cut
sub _vcf_frequency_filtering {
  my $self = shift;
  my $tva = shift;
  my $allele = $tva->base_variation_feature->alt_alleles;  
  my $vf = $tva->base_variation_feature;
  # get the lowest frequency threshold. Threshold can be different for monoallelic and biallelic genes.
  my $frequency_threshold = $self->{config}->{frequency_threshold}; 
  my $vf_cache_name =  $self->{vf_cache_name};
  my $alt_alleles = $allele->[0];
  foreach my $vcf_collection (@{$self->{config}->{vcf_collections}}) {
    my @alleles = grep {$_->allele eq $alt_alleles} @{$vcf_collection->get_all_Alleles_by_VariationFeature($vf)};
    # As soon as we find a frequency which is higher than the frequency_threshold,
    # and variant is not on variant_include_list we can stop.
    my @frequencies = grep {$_->frequency > $frequency_threshold} @alleles;
    if (scalar @frequencies > 0 && !$self->{g2p_vf_cache}->{$vf_cache_name}->{is_on_variant_include_list}) {
      return 0;
    } else {
      $self->_dump_existing_vf_vcf(\@alleles) if (scalar @alleles); 
    }
  }
  return 1;
}

=head2 _dump_existing_vf_vcf

  Arg [1]    : Arrayref $alleles
  Description: Write allele frequencies from VCF files to the log file as G2P_frequencies.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub _dump_existing_vf_vcf {
  my $self = shift;
  my $alleles = shift;
  my @frequencies = map {$_->population->name . '=' . $_->frequency} @$alleles;
  my @sorted_frequencies = sort { $a->frequency <=> $b->frequency } @$alleles;
  $self->highest_frequency($self->{vf_cache_name}, $sorted_frequencies[-1]->frequency);
  $self->write_report('G2P_frequencies', $self->{vf_cache_name}, \@frequencies);
}

=head2 highest_frequency

  Arg [1]    : String $vf_cache_name
  Arg [2]    : Float $frequency
  Description: Getter and setter for highest observed frequency for the current variant in the internal cache which is used
               for filtering later.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub highest_frequency {
  my $self = shift;
  my $vf_cache_name = shift;
  my $f = shift;
  if (defined $vf_cache_name && defined $f) {
    my $highest_frequency = $self->{highest_frequencies}->{$vf_cache_name};
    if (defined $highest_frequency && $highest_frequency < $f || ! defined  $highest_frequency) {
      $self->{highest_frequencies}->{$vf_cache_name} = $f;
    }
  }
  return $self->{highest_frequencies}->{$vf_cache_name};
}

=head2 dump_vf_annotations

  Arg [1]    : TranscriptVariationAllele $tva 
  Description: Write all variation feature related annotations to the log file as G2P_tva_annotations.
               Write flag is_on_variant_include_list to the log file if the variant is on the variant include list. 
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub dump_vf_annotations {
  my $self = shift;
  my $tva = shift;
  my @consequence_types = map { $_->SO_term } @{$tva->get_all_OverlapConsequences};
  my $vf = $tva->base_variation_feature;
  my $allele = $tva->base_variation_feature->alt_alleles;
  my $start = $vf->{start};
  my $end = $vf->{end};

  my $individual = $vf->{individual};
  my $vf_name = $vf->variation_name;
  my $vf_cache_name = $self->{vf_cache_name};
  my $allele_string = $vf->{allele_string};
  my @alleles = split('/', $allele_string);
  my $ref = $alleles[0];
  my $alt = $alleles[1];
  my $seq_region_name = $vf->{chr};

  my $is_on_variant_include_list = $self->{g2p_vf_cache}->{$vf_cache_name}->{is_on_variant_include_list} || 0;

  my $params = $self->{user_params};
  my $tr = $tva->transcript;
  my $refseq = $tr->{_refseq} || 'NA';
  my $hgvs_t = $tva->hgvs_transcript || 'NA';
  my $hgvs_p = $tva->hgvs_protein || 'NA';

  my $pph_score   = (defined $tva->polyphen_score) ? $tva->polyphen_score : 'NA';
  my $pph_pred    = (defined $tva->polyphen_prediction) ? $tva->polyphen_prediction : 'NA';
  my $sift_score  = (defined $tva->sift_score) ? $tva->sift_score : 'NA';
  my $sift_pred   = (defined $tva->sift_prediction) ? $tva->sift_prediction : 'NA';

  my $g2p_data = {
    'vf_name' => $vf_name,
    'is_on_variant_include_list' => $is_on_variant_include_list,
    'transcript_stable_id' => $tr->stable_id,
    'consequence_types' => join(',', @consequence_types),
    'refseq' => $refseq,
    'hgvs_t' => $hgvs_t,
    'hgvs_p' => $hgvs_p,
    'vf_location' => "$seq_region_name:$start-$end $ref/$alt",
    'sift_score' => "$sift_score",
    'sift_prediction' => $sift_pred,
    'polyphen_score' => "$pph_score",
    'polyphen_prediction' => $pph_pred,
  };
  $self->write_report('G2P_tva_annotations', $vf_cache_name, $tr->stable_id, $g2p_data);
  $self->write_report('is_on_variant_include_list', $vf_cache_name) if ($is_on_variant_include_list);
}

=head2 dump_individual_annotations

  Arg [1]    : TranscriptVariationAllele $tva
  Description: Write individual specific information to the log file as G2P_individual_annotations.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub dump_individual_annotations {
  my $self = shift;
  my $tva = shift;
  my $zyg = shift;
  my $vf = $tva->base_variation_feature;
  my $individual = $vf->{individual};
  my $vf_cache_name = $self->{vf_cache_name};
  my $transcript = $tva->transcript;
  my $transcript_stable_id = $transcript->stable_id;
  my $gene_stable_id = get_gene_stable_id($transcript);
  $self->write_report('G2P_individual_annotations', join("\t", $gene_stable_id, $transcript_stable_id, $vf_cache_name, $zyg, $individual));
}

=head2 read_gene_data_from_file

  Arg [1]    : String $file
  Description: Read panel data from file into the internal cache. Extract gene symbol, gene symbol synonyms or previously assigned symbols,
               allelic requirement and gene-disease confidence values.
               Get G2P CSV dump from https://www.ebi.ac.uk/gene2phenotype/downloads.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub read_gene_data_from_file {
  my $self = shift;
  my $file = shift;
  my (@headers, %gene_data);

  my $assembly =  $self->{config}->{assembly};
  die("ERROR: No file specified or could not read from file ".($file || '')."\n") unless $file && -e $file;

  my @confidence_levels = @{$self->{user_params}->{confidence_levels}};

  # determine file type
  my $file_type;
  my $fh = FileHandle->new($file, 'r');
  while (<$fh>) {
    chomp;
      if (/Model_Of_Inheritance/) {
        $file_type = 'panelapp';
      } elsif (/allelic requirement/) {
        $file_type = 'g2p';
      } else {
        $file_type = 'unknown';
      }
      last;
  }
  $fh->close();
  if ($file_type eq 'unknown') {
    if ($file =~ /gz$/) { 
      die("ERROR: G2P plugin can only read uncompressed data\n");
    } else {
      die("ERROR: Could not recognize input file format. Format must be one of panelapp, g2p or custom. Check website for details: https://www.ebi.ac.uk/gene2phenotype/g2p_vep_plugin\n");
    }
  }

  if ($file_type eq 'panelapp') {
    my @headers = ();
    my $csv = Text::CSV->new ({ sep_char => "\t" });
    open my $fh, "<:encoding(utf8)", "$file" or die "$file: $!";
    while ( my $row = $csv->getline( $fh ) ) {
      unless (@headers) {
        @headers = @$row;
      } else {
        my %tmp = map {$headers[$_] => $row->[$_]} (0..$#headers);
        my $gene_symbol = $tmp{"Gene Entity Symbol"};
        my $ensembl_gene_id = "";
        if ($assembly eq 'GRCh37') { 
          $ensembl_gene_id = $tmp{"EnsemblId(GRch37)"};
        } else { # GRCh38
          $ensembl_gene_id = $tmp{"EnsemblId(GRch38)"};
        }
        if ($ensembl_gene_id) {
          foreach my $column_name ('Gene Entity Symbol', 'Gene Symbol') {
            if (defined $tmp{$column_name}) {
              push @{$gene_data{$ensembl_gene_id}->{"gene_xrefs"}}, $tmp{$column_name};
              $self->write_report('G2P_list', $tmp{$column_name});
            }
          }
          my @ars = ();
          my $allelic_requirement_panel_app = $tmp{"Model_Of_Inheritance"};
          if ($allelic_requirement_panel_app =~ m/BOTH monoallelic and biallelic/) {
            push @ars, 'monoallelic';
            push @ars, 'biallelic';
          } elsif ($allelic_requirement_panel_app eq 'X-LINKED: hemizygous mutation in males, biallelic mutations in females') {
            push @ars, 'hemizygous';
          } elsif ($allelic_requirement_panel_app eq 'X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)') {
            push @ars, 'x-linked dominant';
          } elsif ($allelic_requirement_panel_app =~ m/MONOALLELIC/) {
            push @ars, 'monoallelic';
          } elsif ($allelic_requirement_panel_app =~ m/BIALLELIC/) {
            push @ars, 'biallelic';
          } elsif ($allelic_requirement_panel_app =~ m/MITOCHONDRIAL/) {
            push @ars, 'mitochondrial';
          } else {
            $self->write_report('log', "no allelelic_requirement for $ensembl_gene_id");
          }
          foreach my $ar (@ars) {
            push @{$gene_data{$ensembl_gene_id}->{"allelic requirement"}}, $ar;
          }
        } else {
          $self->write_report('log', "no ensembl gene id");
        }
      }
    }
    $csv->eof or $csv->error_diag();
    close $fh;
  }

  if ($file_type eq 'g2p') {
    # this regexp allows for nested ",", e.g.
    # item,description
    # cheese,"salty,delicious"
    my $re = qr/(?: "\( ( [^()""]* ) \)" |  \( ( [^()]* ) \) |  " ( [^"]* ) " |  ( [^,]* ) ) , \s* /x;

    my $fh = FileHandle->new($file, 'r');

    while(<$fh>) {
      chomp;
      $_ =~ s/\R//g;
      my @split = grep defined, "$_," =~ /$re/g;
      unless(@headers) {
        if ($file_type eq 'g2p') {
          @headers = map {s/\"//g; $_} @split;
        } else {
          @headers = @split;
        }
      }
      else {
        my %tmp = map {$headers[$_] => $split[$_]} (0..$#split);
        die("ERROR: Gene symbol column not found\n$_\n") unless $tmp{"gene symbol"};
        $self->write_report('G2P_list', $tmp{"gene symbol"}, $tmp{"hgnc id"}, $tmp{"DDD category"});
        my $confidence_category = $tmp{"DDD category"} || $tmp{"confidence category"}; # deprecate use of DDD category
        next if (!grep{$_ eq $confidence_category} @confidence_levels);
        my $gene_symbol = $tmp{"gene symbol"};
        push @{$gene_data{$gene_symbol}->{"gene_xrefs"}}, split(';', $tmp{"prev symbols"}) if defined($tmp{"prev symbols"});
        push @{$gene_data{$gene_symbol}->{"gene_xrefs"}}, split(';', $tmp{"previous gene symbols"}) if defined($tmp{"previous gene symbols"});
        push @{$gene_data{$gene_symbol}->{"gene_xrefs"}}, $tmp{"gene symbol"};
        push @{$gene_data{$gene_symbol}->{"HGNC"}}, $tmp{"hgnc id"};
        push @{$gene_data{$gene_symbol}->{"confidence_category"}}, $confidence_category;
        push @{$gene_data{$gene_symbol}->{"allelic requirement"}}, $tmp{"allelic requirement"} if ($tmp{"allelic requirement"});
        push @{$gene_data{$gene_symbol}->{'confidence_value'}}, $tmp{"confidence value flag"}  if ($tmp{"confidence value flag"}); # adding the confidence value flag "Requires clinical review if defined"
      }
    }
    $fh->close;
  }
  return \%gene_data;
}

=head2 gene_data

  Arg [1]    : String $gene_symbol
  Example    : $gene_data = $self->gene_data('PRKAR1A')
  Description: Get all panel specific data for the given gene symbol.
               TODO add example of gene_data hash
  Returntype : Hashref $gene_data, for example:
                  {
                    'gene_xrefs' => [
                      'PEPD'
                    ],
                    'allelic requirement' => [
                      'biallelic'
                    ]
                  };
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub gene_data {
  my ($self, $gene_symbol) = @_;
  my $gene_data = $self->{gene_data}->{$gene_symbol};
  if (!$gene_data) {
    my $prev_gene_symbol = $self->{prev_symbol_mappings}->{$gene_symbol};
    return $prev_gene_symbol ? $self->{gene_data}->{$prev_gene_symbol} : undef;
  } 
  return $gene_data;
}


=head2 synonym_mappings

  Description: Create a hashref in the internal cache under the key prev_symbol_mappings which
               maps any previously assigned gene symbols to the latest gene symbol.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub synonym_mappings {
  my $self = shift;
  my $gene_data = $self->{gene_data};
  my $synonym_mappings = {};
  foreach my $gene_symbol (keys %$gene_data) {
    foreach my $prev_symbol (@{$gene_data->{$gene_symbol}->{'gene_xrefs'}}) {
      $synonym_mappings->{$prev_symbol} = $gene_symbol;
    }
  }
  $self->{prev_symbol_mappings} = $synonym_mappings;
}

sub hgnc_mappings {
  my $self = shift;
  my $gene_data = $self->{gene_data};
  my $hgnc_mapping = {};
  foreach my $gene_symbol (keys %$gene_data){
    foreach my $hgnc_id (@{$gene_data->{$gene_symbol}->{'HGNC'}}) {
      $hgnc_mapping->{$hgnc_id} = $gene_symbol;
    }
  }
  $self->{hgnc_mapping} = $hgnc_mapping;
  
}

=head2 write_report

  Arg [1]    : String $flag
  Arg [2]    : Array of values that need to be written to the log file together with the flag.
  Description: write_report is called to write results and annotations to the log file.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub write_report {
  my $self = shift;
  my $flag = shift;
  my $log_dir = $self->{user_params}->{log_dir};
  my $log_file = "$log_dir/$$.txt";
  open(my $fh, '>>', $log_file) or die "Could not open file '$flag $log_file' $!\n";
  if ($flag eq 'G2P_list') {
    my ($gene_symbol, $DDD_category) = @_;
    $DDD_category ||= 'Not assigned';
    print $fh "$flag\t$gene_symbol\t$DDD_category\n";
  } elsif ($flag eq 'G2P_in_vcf') {
    my $gene_symbol = shift;
    print $fh "$flag\t$gene_symbol\n";
  } elsif ($flag eq 'G2P_complete') {
    print $fh join("\t", $flag, @_), "\n";
  } elsif ($flag eq 'log') {
    print $fh join("\t", $flag, @_), "\n";
  } elsif ($flag eq 'is_on_variant_include_list') {
    my ($vf_name) = @_;
    print $fh "$flag\t$vf_name\n";
  } elsif ($flag eq 'G2P_gene_data') {
    my ($gene_id, $gene_data, $gene_xrefs, $hgnc, $con_category, $con_flag) = @_;
    my $ar = join(',', @{$gene_data->{'allelic requirement'}});
    my %seen;
    $seen{$_} = 1 foreach @{$gene_xrefs};
    my @unique = keys %seen;
    my $xrefs = join(',', grep {$_ !~ /^ENS/} sort @unique);
    my $hgnc_id = "HGNC:".@{$hgnc}[0].";" if (defined $hgnc);
    my $conf = "Note:".@{$con_flag}[0] if (defined $con_flag);
    my $category = "G2P confidence category:".@{$con_category}[0] if defined ($con_category);
    # if statements to determine what should be printed based on $hgnc and $conf 
    print $fh join("\t", $flag, $gene_id, $ar, $xrefs, $hgnc_id, $category, $conf), "\n"  if (defined $hgnc_id && defined $conf && defined $category);
    print $fh join("\t", $flag, $gene_id, $ar, $xrefs, $category, $conf), "\n" if (!defined $hgnc_id && defined $conf && defined $category);
    print $fh join("\t", $flag, $gene_id, $ar, $xrefs, $hgnc_id, $conf), "\n" if (defined $hgnc_id && defined $conf);
    print $fh join("\t", $flag, $gene_id, $ar, $xrefs, $hgnc_id, $category), "\n"  if (defined $hgnc_id && !defined $conf & defined $category);
    print $fh join("\t", $flag, $gene_id, $ar, $xrefs, $category), "\n"  if (!defined $hgnc_id && !defined $conf && defined $category);
    print $fh join("\t", $flag, $gene_id, $ar, $xrefs), "\n"  if (!defined $hgnc_id && !defined $conf && !defined $category);
  } elsif ($flag eq 'G2P_frequencies') {
    my ($vf_name, $frequencies) = @_;
    print $fh join("\t", $flag, $vf_name, join(',', @$frequencies)), "\n";
  } elsif ($flag eq 'G2P_tva_annotations') {
    my ($vf_name, $transcript_stable_id, $data) = @_;
    $data = join(';', map {"$_=$data->{$_}"} sort keys %$data);
    print $fh join("\t", $flag, $vf_name, $transcript_stable_id, $data), "\n";
  } elsif ($flag eq 'G2P_existing_vf_annotations') {
    my ($vf_name, $data) = @_;
    $data = join(';', map {"$_=$data->{$_}"} sort keys %$data);
    print $fh join("\t", $flag, $vf_name, $data), "\n";
  } elsif ($flag eq 'G2P_transcript_data') {
    print $fh join("\t", $flag, @_), "\n";
  } elsif ($flag eq 'G2P_individual_annotations') { 
    print $fh join("\t", $flag, @_), "\n";
  } else {
    print $fh "Did not recognize flag, @_\n";
  }
  close $fh;
}

=head2 finish

  Description: Call generate_report subroutine.
  Returntype : None
  Exceptions : None
  Caller     : Called by VEP::Runner after it completed the annotation of all variants from the input VCF file.
  Status     : Stable

=cut
sub finish {
  my $self = shift;
  $self->generate_report;
}

=head2 generate_report

  Description: Parses the log files into a hashref $result_summary.
               Calls html_and_txt_data($result_summary) to create specific hashrefs for generating the TXT and HTML output files.
  Returntype : None
  Exceptions : None
  Caller     : G2P::finish
  Status     : Stable

=cut
sub generate_report {
  my $self = shift;
  my $result_summary = $self->parse_log_files;
  my ($html_data, $txt_data) = $self->html_and_txt_data($result_summary);
  $self->write_html_output($result_summary, $html_data, $result_summary->{canonical_transcripts}, $result_summary->{gene_xrefs});
  $self->write_txt_output($txt_data, $result_summary->{gene_xrefs});
}

=head2 write_txt_output

  Arg [1]    : Hashref $txt_output_data
  Arg [2]    : Hashref $gene_xrefs
  Description: Genereates the TXT output file.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub write_txt_output {
  my $self = shift;
  my $txt_output_data = shift; 
  my $gene_xrefs = shift;
  my $txt_output_file = $self->{user_params}->{txt_report};
  my $fh_txt = FileHandle->new($txt_output_file, 'w');
  foreach my $individual (sort keys %$txt_output_data) {
    foreach my $gene_id (keys %{$txt_output_data->{$individual}}) {
      my $gene_id_title = (defined $gene_xrefs->{$gene_id}) ? "$gene_id(" .  $gene_xrefs->{$gene_id} . ")" : $gene_id;
      # Loop over all observed allelic requirements which are fulfilled based on the number
      # of variants that have been found in the input VCF
      foreach my $ar (keys %{$txt_output_data->{$individual}->{$gene_id}}) {
        foreach my $tr_stable_id (keys %{$txt_output_data->{$individual}->{$gene_id}->{$ar}}) {
          my $is_canonical = $txt_output_data->{$individual}->{$gene_id}->{$ar}->{$tr_stable_id}->{is_canonical};
          my $canonical_tag = ($is_canonical) ? 'is_canonical' : 'not_canonical';
          my $req =  $txt_output_data->{$individual}->{$gene_id}->{$ar}->{$tr_stable_id}->{REQ};
          my $variants = join(';', @{$txt_output_data->{$individual}->{$gene_id}->{$ar}->{$tr_stable_id}->{variants}});
          print $fh_txt join("\t", $individual, $gene_id_title, $tr_stable_id, $canonical_tag, "OBS=$ar", "REQ=$req", $variants), "\n";
        }
      }
    }
  }
  $fh_txt->close();
}

=head2 write_html_output

  Arg [1]    : Hashref $result_summary
  Arg [2]    : Hashref $html_data
  Arg [3]    : Hashref $canonical_transcripts
  Arg [4]    : Hashref $gene_xrefs
  Description: Generates the HTML output file.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub write_html_output {
  my $self = shift;
  my $result_summary = shift;
  my $html_data = shift;
  my $canonical_transcripts = shift;
  my $gene_xrefs = shift;

# G2P genes in G2P input panel file
  my $count_g2p_genes = keys %{$result_summary->{g2p_list}};
# G2P genes in input VCF file: How many G2P genes overlap any of the variants from the input VCF.
  my $count_in_vcf_file = keys %{$result_summary->{in_vcf_file}};
# G2P complete genes in input VCF file
  my $count_g2p_complete_genes = scalar keys %{$result_summary->{g2p_complete_genes}};

  my @frequencies_header = (); 

  foreach my $short_name (sort keys %{$self->{population_names}}) {
    my $text = $af_key_2_population_name->{$short_name} || 'No description';
    push @frequencies_header, "<a style=\"cursor: pointer\" data-placement=\"top\" data-toggle=\"tooltip\" data-container=\"body\" title=\"$text\">$short_name</a>";
  }

  my $count = 1;
  my @new_header = (
    'Variant location and alleles (REF/ALT)',
    'Variant name (* indicates that variant is on variant include list)',
    'Existing name', 
    'Zygosity', 
    'All allelic requirements from G2P DB',
    'Consequence types', 
    'ClinVar annotation', 
    'SIFT', 
    'PolyPhen', 
    'Novel variant', 
    'Has been failed by Ensembl', 
    @frequencies_header,
    'HGVS transcript', 
    'HGVS protein', 
    'RefSeq IDs', 
  );

  my $html_output_file = $self->{user_params}->{html_report};
  my $fh_out = FileHandle->new($html_output_file, 'w');
  print $fh_out stats_html_head();
  print $fh_out "<div class='main_content container'>";

 
  print $fh_out "<h1>G2P report</h1>";
  print $fh_out "<p>Input and output files:</p>";

  print $fh_out "<dl class='dl-horizontal'>";
  print $fh_out "<dt>G2P list</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{file} .  "</dd>";
  print $fh_out "<dt>Log directory</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{log_dir} .  "</dd>";
  print $fh_out "<dt>HTML report</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{html_report} .  "</dd>";
  print $fh_out "<dt>TXT report</dt>";
  print $fh_out "<dd>" . $self->{user_params}->{txt_report} .  "</dd>";
  print $fh_out "</dl>";

  print $fh_out "<p>Counts:</p>";
  print $fh_out "<dl class='dl-horizontal text-overflow'>";
  print $fh_out "<dt>$count_g2p_genes</dt>";
  print $fh_out "<dd>G2P genes</dd>";
  print $fh_out "<dt>$count_in_vcf_file</dt>";
  print $fh_out "<dd>G2P genes in input VCF file</dd>";
  print $fh_out "<dt>$count_g2p_complete_genes</dt>";
  print $fh_out "<dd>G2P complete genes in input VCF file</dd>";
  print $fh_out "</dl>";

  
  print $fh_out "<h1>Summary of G2P complete genes per individual</h1>";
  print $fh_out "<p>G2P complete gene: A sufficient number of variant hits for the observed allelic requirement in at least one of the gene's transcripts. Variants are filtered by frequency.</p>";
  print $fh_out "<p>G2P gene was filtered using HGNC_id</p>" if (!defined $self->{user_params}->{filter_by_gene_symbol});
  print $fh_out "<p>Frequency thresholds and number of required variant hits for each allelic requirement:</p>";

  print $fh_out "<table class='table table-bordered'>";
  print $fh_out "<thead>";
  print $fh_out "<tr><th>Allelic requirement</th><th>Frequency threshold for filtering</th><th>Variant counts by zygosity</th></tr>";
  print $fh_out "</thead>";
  print $fh_out "<tbody>";
  foreach my $ar (sort keys %$allelic_requirements) {
    my $af = $allelic_requirements->{$ar}->{af};
    my $rules =  $allelic_requirements->{$ar}->{rules};
    my $rule = join(' OR ', map {"$_ >= $rules->{$_}"} keys %$rules);
    print $fh_out "<tr><td>$ar</td><td>$af</td><td>$rule</td></tr>";
  }
  print $fh_out "</tbody>";
  print $fh_out "</table>";


my $switch =<<SHTML;
<form>
<div class="checkbox">
  <label>
    <input class="target" type="checkbox"> Show only canonical transcript
  </label>
</div>
</form>
SHTML

  print $fh_out $switch if (!$self->{user_params}->{only_mane}); 
 
  foreach my $individual (sort keys %$html_data) {
    foreach my $gene_id (keys %{$html_data->{$individual}}) {
      my $gene_id_title = (defined $gene_xrefs->{$gene_id}) ? "$gene_id(" .  $gene_xrefs->{$gene_id} . ")" : $gene_id;
      foreach my $ar (keys %{$html_data->{$individual}->{$gene_id}}) {
        print $fh_out "<ul>\n";
        foreach my $transcript_stable_id (keys %{$html_data->{$individual}->{$gene_id}->{$ar}}) {
          my $class = ($canonical_transcripts->{$transcript_stable_id}) ? 'is_canonical' : 'not_canonical';
          print $fh_out "<li><a class=\"$class\" href=\"#$individual\_$gene_id_title\_$ar\_$transcript_stable_id\">" . "$individual &gt; $gene_id_title &gt; $ar &gt; $transcript_stable_id" . "</a> </li>\n";
        }
        print $fh_out "</ul>\n";
      }
    }
  }

  foreach my $individual (sort keys %$html_data) {
    foreach my $gene_id (keys %{$html_data->{$individual}}) {
      my $gene_id_title = (defined $gene_xrefs->{$gene_id}) ? "$gene_id(" .  $gene_xrefs->{$gene_id} . ")" : $gene_id;
      foreach my $ar (keys %{$html_data->{$individual}->{$gene_id}}) {
        foreach my $transcript_stable_id (keys %{$html_data->{$individual}->{$gene_id}->{$ar}}) {
          my $class = ($canonical_transcripts->{$transcript_stable_id}) ? 'is_canonical' : 'not_canonical';
          print $fh_out "<div class=\"$class\">";
          my $name = "$individual\_$gene_id_title\_$ar\_$transcript_stable_id";
          my $title = "$individual &gt; $gene_id_title &gt; $ar &gt; $transcript_stable_id";
          print $fh_out "<h3><a name=\"$name\"></a>$title <a title=\"Back to Top\" data-toggle=\"tooltip\" href='#top'><span class=\"glyphicon glyphicon-arrow-up\" aria-hidden=\"true\"></span></a></h3>\n";
          print $fh_out "<div class=\"table-responsive\" style=\"width:100%\">\n";
          print $fh_out "<TABLE  class=\"table table-bordered table-condensed\" style=\"margin-left: 2em\">";
          print $fh_out "<thead>\n";
          print $fh_out "<tr>" . join('', map {"<th>$_</th>"} @new_header) . "</tr>\n";
          print $fh_out "</thead>\n";
          print $fh_out "<tbody>\n";
          foreach my $vf_data (@{$html_data->{$individual}->{$gene_id}->{$ar}->{$transcript_stable_id}}) {
            my $data_row = $vf_data->[0];
            my @tds = ();
            foreach my $cell (@$data_row) {
              my $value = $cell->[0];
              my $class = $cell->[1];
              if ($class) {
                push @tds, "<td class=\"$class\">$value</td>";
              } else {
                push @tds, defined $value ? "<td>$value</td>" : "<td>NA</td>";
              }
            }
            print $fh_out "<tr>", join('', @tds), "</tr>\n";
          }
          print $fh_out "</tbody>\n";
          print $fh_out "</TABLE>\n";
          print $fh_out "</div>\n";
          print $fh_out "</div>\n";
        }
      }
    }
  }
  print $fh_out stats_html_tail();
}

=head2 html_and_txt_data

  Arg [1]    : Hashref $result_summary
  Description: Create two hashrefs  $html_output_data and $txt_output_data which will be used to create the HTML and TXT output files.
  Returntype : Array of Hashref $html_output_data, Hashref $txt_output_data
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub html_and_txt_data {
  my $self = shift;
  my $result_summary = shift;
  my $results_by_individual = $result_summary->{results_by_individual};

  my $tva_annotation_data = $result_summary->{tva_annotation_data};
  my $vf_annotation_data = $result_summary->{vf_annotation_data};
  my $frequency_data = $result_summary->{frequency_data};
  my $canonical_transcripts = $result_summary->{canonical_transcripts};
  my $gene2ar = $result_summary->{gene2ar}; # All allelic requirements that have been reported for the  gene in the G2P database

  my @frequencies_header = sort keys %{$self->{population_names}};

  my $assembly = $self->{config}->{assembly};
  my $html_output_data = {};
  my $txt_output_data = {};

  my $prediction2bgcolor = {
    'probably damaging' => 'danger',
    'deleterious' => 'danger',
    'possibly damaging' => 'warning',
    'unknown'  => 'warning',
    'benign' => 'success',
    'tolerated' => 'success',
  };

  foreach my $individual (sort keys %$results_by_individual) {

    foreach my $gene_id (keys %{$results_by_individual->{$individual}}) {
      # All required allelic requirements as reported in the G2P database
      my $required_allelic_requirement = join(',', keys %{$gene2ar->{$gene_id}});
      # Loop over all allelic requirements that are fulfilled based on the number of
      # variants that have been found.
      foreach my $ar (keys %{$results_by_individual->{$individual}->{$gene_id}}) {
        foreach my $transcript_stable_id (keys %{$results_by_individual->{$individual}->{$gene_id}->{$ar}}) {
          my $zyg2vf = $results_by_individual->{$individual}->{$gene_id}->{$ar}->{$transcript_stable_id}; 
          foreach my $zygosity (keys %$zyg2vf) {
            foreach my $vf_name (@{$zyg2vf->{$zygosity}}) {
              my $tva_data = $tva_annotation_data->{$vf_name}->{$transcript_stable_id};
              my $vf_data = $vf_annotation_data->{$vf_name}; 
              if (!$vf_data) {
                print STDERR "No vf_data for: $vf_name\n"; 
              } 
              my $hash = {};
              foreach my $pair (split /;/, defined $vf_data ? "$tva_data;$vf_data" : $tva_data) {
                my ($key, $value) = split('=', $pair, 2);
                $value ||= '';
                $hash->{$key} = $value;
              }
              my $vf_location = $hash->{vf_location};
              my $existing_name = $hash->{existing_name};
              my $ensembl_url;
              if($assembly eq 'GRCh37') {
                $ensembl_url = "https://grch37.ensembl.org";
              }
              else {
                $ensembl_url = "https://ensembl.org";
              }
              if (defined($existing_name) && $existing_name ne 'NA' ) {
                $existing_name = "<a href=\"$ensembl_url/Homo_sapiens/Variation/Explore?v=$existing_name\">$existing_name</a>";
              }
              my $is_on_variant_include_list = $hash->{is_on_variant_include_list};
              my $refseq = $hash->{refseq};
              my $failed = $hash->{failed};
              my $clin_sign = $hash->{clin_sig};
              my $novel = $hash->{novel};
              my $hgvs_t = $hash->{hgvs_t};
              my $hgvs_p = $hash->{hgvs_p};
              my $consequence_types = $hash->{consequence_types};
              my $sift_score = $hash->{sift_score} || '0.0';
              my $sift_prediction = $hash->{sift_prediction};
              my $sift = 'NA';
              my $sift_class = '';
              if ($sift_prediction ne 'NA') {
                $sift = "$sift_prediction(" . "$sift_score)";
                $sift_class = $prediction2bgcolor->{$sift_prediction};
              }
              my $polyphen_score = $hash->{polyphen_score} || '0.0';
              my $polyphen_prediction = $hash->{polyphen_prediction};
              my $polyphen = 'NA';
              my $polyphen_class = '';
              if ($polyphen_prediction ne 'NA') {
                $polyphen = "$polyphen_prediction($polyphen_score)";
                $polyphen_class =  $prediction2bgcolor->{$polyphen_prediction};
              }
              
              $hash->{frequencies} = join(',', keys %{$frequency_data->{$vf_name}}) || 'NA';
              my %frequencies_hash = ();
              if ($hash->{frequencies} ne 'NA') {
                %frequencies_hash = split /[,=]/, $hash->{frequencies};
              }
              my @frequencies = ();
              my @txt_output_frequencies = ();
              foreach my $population (@frequencies_header) {
                my $frequency = $frequencies_hash{$population} || '';
                push @frequencies, ["$frequency"];
                if ($frequency) {
                  push @txt_output_frequencies, "$population=$frequency";
                }
              }
              my $is_canonical = ($canonical_transcripts->{$transcript_stable_id}) ? 1 : 0;
              my ($location, $alleles) = split(' ', $vf_location);
              $location =~ s/\-/:/;
              $alleles =~ s/\//:/;
              $vf_name .= "*" if ($is_on_variant_include_list);
              push @{$html_output_data->{$individual}->{$gene_id}->{$ar}->{$transcript_stable_id}}, [[
                [$vf_location], 
                [$vf_name], 
                [$existing_name], 
                [$zygosity], 
                [$required_allelic_requirement],
                [$consequence_types], 
                [$clin_sign], 
                [$sift, $sift_class], 
                [$polyphen, $polyphen_class], 
                [$novel], 
                [$failed], 
                @frequencies,
                [$hgvs_t], 
                [$hgvs_p], 
                [$refseq] 
              ], $is_canonical];

              my $txt_output_variant = "$location:$alleles:$zygosity:$consequence_types:SIFT=$sift:PolyPhen=$polyphen";
              if (@txt_output_frequencies) {
                $txt_output_variant .= ':' . join(',', @txt_output_frequencies);
              }
              $txt_output_data->{$individual}->{$gene_id}->{$ar}->{$transcript_stable_id}->{is_canonical} = $is_canonical;
              $txt_output_data->{$individual}->{$gene_id}->{$ar}->{$transcript_stable_id}->{REQ} = $required_allelic_requirement;
              push @{$txt_output_data->{$individual}->{$gene_id}->{$ar}->{$transcript_stable_id}->{variants}}, $txt_output_variant;
            }
          }
        }
      }
    }
  }
  return ($html_output_data, $txt_output_data);
}

=head2 parse_log_files

  Description: Read content from log files and write data into hashrefs which will
               be used to extract G2P complete genes and write results into TXT and
               HTML output files.
  Returntype : Hashref of hashrefs
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub parse_log_files {
  my $self = shift;

  my $log_dir = $self->{user_params}->{log_dir}; 
  my @files = <$log_dir/*>;
  my $individual_data = {};
  my $frequency_data = {};
  my $vf_annotation_data = {};
  my $tva_annotation_data = {};
  my $canonical_transcripts = {};
  my $all_g2p_genes = {};
  my $vcf_g2p_genes = {};
  my $ar_data = {}; # all allelic requirements that have been reported for the gene in the G2P database
  my $g2p_transcripts = {};
  my $gene_xrefs = {};
  

  foreach my $file (@files) {
    my $fh = FileHandle->new($file, 'r');
    while (<$fh>) {
     chomp;
      next if /^log/;
      if (/^G2P_list/) {
        my ($flag, $gene_id, $DDD_category) = split/\t/;
        $all_g2p_genes->{$gene_id} = 1;
      }
      #G2P_individual_annotations  ENSG00000091140 ENST00000450038 7_107545113_T/C HOM P10
      elsif (/^G2P_individual_annotations/) {
        my ($flag, $gene_stable_id, $transcript_stable_id, $vf_cache_name, $zyg, $individual) = split/\t/;
        $individual_data->{$individual}->{$gene_stable_id}->{$transcript_stable_id}->{$zyg}->{$vf_cache_name} = 1;
      }
      #G2P_frequencies     17_82929274_A/G AA=0.0007289,AFR=0,AMR=0.0014,EA=0.005595,EAS=0,EUR=0.006,SAS=0,gnomAD=0.002923,gnomAD_AFR=0.0005175,gnomAD_AMR=0.002206,gnomAD_ASJ=0.003392,gnomAD_EAS=0,gnomAD_FIN=0.003583,gnomAD_NFE=0.004446,gnomAD_OTH=0.003814,gnomAD_SAS=0.0002618
      elsif (/^G2P_frequencies/) {
        my ($flag, $vf_cache_name, $frequencies) = split/\t/;
        $frequency_data->{$vf_cache_name}->{$frequencies} = 1;
        $self->store_population_names($frequencies);
        my $highest_frequency = get_highest_frequency($frequencies);
        $self->highest_frequency($vf_cache_name, $highest_frequency);
      }
      #G2P_tva_annotations 17_82929274_A/G ENST00000355528 consequence_types=splice_region_variant,intron_variant;hgvs_p=NA;hgvs_t=ENST00000355528.9:c.2852+3A>G;is_on_variant_include_list=0;polyphen_prediction=NA;polyphen_score=NA;refseq=NM_005993.5;sift_prediction=NA;sift_score=NA;transcript_stable_id=ENST00000355528;vf_location=17:82929274-82929274 A/G;vf_name=id_17_82929274_A_G
      elsif (/^G2P_tva_annotations/) {
        my ($flag, $vf_cache_name, $transcript_stable_id, $annotations) = split/\t/;
        $tva_annotation_data->{$vf_cache_name}->{$transcript_stable_id} = $annotations;
      }
      #G2P_existing_vf_annotations 17_82941481_C/T clin_sig=NA;existing_name=rs780053410;failed=no;novel=no
      elsif (/^G2P_existing_vf_annotations/) {
        my ($flag, $vf_cache_name, $annotations) = split/\t/;
        $vf_annotation_data->{$vf_cache_name} = $annotations;
      }
      #G2P_gene_data ENSG00000141556 biallelic TBCD HGNC:15432 if not filter_by_gene_symbol
      elsif (/^G2P_gene_data/) {
        my ($flag, $gene_id, $ars, $xrefs, $hgnc, $con_flag) = split/\t/;
        foreach my $ar (split(',', $ars)) {
          $ar_data->{$gene_id}->{$ar} = 1;
        }
        
        if (defined $hgnc && !defined $self->{user_params}->{filter_by_gene_symbol}){
          my $x_hgnc = $xrefs . " $hgnc  $con_flag" if (defined $con_flag);
          $x_hgnc = $xrefs . " $hgnc " if (!defined $con_flag);
          $gene_xrefs->{$gene_id} = $x_hgnc;
        }
        if (!defined $hgnc && defined $self->{user_params}->{filter_by_gene_symbol} ) {
          $gene_xrefs->{$gene_id} = $xrefs . " $con_flag" if (defined $con_flag) ;
          $gene_xrefs->{$gene_id} = $xrefs if (!defined $con_flag);
        }
     }
      #G2P_in_vcf  ENSG00000141556
      elsif (/^G2P_in_vcf/) {
        my ($flag, $gene_id) = split/\t/;
        $vcf_g2p_genes->{$gene_id} = 1;
      }
      #G2P_transcript_data ENSG00000141556 ENST00000355528 is_canonical
      elsif (/^G2P_transcript_data/) {
        my ($flag, $gene_id, $transcript_id, $is_canonical) = split/\t/;
        $canonical_transcripts->{$transcript_id} = 1;
      }
      elsif (/^is_on_variant_include_list/) {
        my ($flag, $vf_cache_name) =  split/\t/;
        $self->{g2p_vf_cache}->{$vf_cache_name}->{is_on_variant_include_list} = 1;
      }
    }
    $fh->close;
  }
  
  my $results_by_individual = $self->get_results_by_individual($individual_data, $ar_data);
  my $g2p_complete_genes = get_g2p_complete_genes($results_by_individual);

  return {
    frequency_data => $frequency_data,
    vf_annotation_data => $vf_annotation_data,
    tva_annotation_data => $tva_annotation_data,
    canonical_transcripts => $canonical_transcripts,
    results_by_individual => $results_by_individual,
    gene2ar => $ar_data,
    gene_xrefs => $gene_xrefs,
    in_vcf_file => $vcf_g2p_genes,
    g2p_complete_genes => $g2p_complete_genes,
    g2p_list => $all_g2p_genes,
  };
}

=head2 get_g2p_complete_genes

  Arg [1]    : Hashref $results_by_individual
  Example    : my $g2p_complete_genes = get_g2p_complete_genes($results_by_individual);
  Description: Extract all gene ids from the $results_by_individual hashref.
  Returntype : Hashref $g2p_complete_genes
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub get_g2p_complete_genes {
  my $results_by_individual = shift;
  my $g2p_complete_genes = {};
  foreach my $individual (keys %$results_by_individual) {
    foreach my $gene_id (keys %{$results_by_individual->{$individual}}) {
      $g2p_complete_genes->{$gene_id} = 1;
    }
  }
  return $g2p_complete_genes;
}

=head2 get_results_by_individual

  Arg [1]    : Hashref $individual_data
  Arg [2]    : Hashref $ar_data
  Example    : $results_by_individual = $self->get_results_by_individual(
               {
                'person_a' => {                # individual id
                  'PRDM15' => {                # gene id
                    'ENST00000398548' => {     # transcript id
                      'HOM' => {               # zygosity: HET, HOM
                        '21_41839757_C/T' => 1 # variant name
                      }
                    }
                  }
                }
               }, {
                'PRDM15' => {
                  'monoallelic' => 1,
                },
               });  
  Description: Creates a new hashref $results_by_individual which orders data from the log file by individual and gene.
               For each gene it then list all transcripts and variants that pass filtering under the respective allelic
               requirement. 
  Returntype : Hashref $results_by_individual
               Here is an example, which is made up. For the monoallelic gene PRDM15 in person a we found a homozygous variant
               in transcript ENST00000398548 which passes all filters and is reported in the result set.
               {
                 'person_a' => {              # individual id 
                   'PRDM15' => {              # gene id
                     'monoallelic' => {       # allelic requirement, contains all variants that fulfil allelic requirement
                       'ENST00000398548' => { # transcript id
                         'HOM' => [           # zygosity
                           '21_41839757_C/T'  # variants that passed filtering by frequency threshold
                         ]
                       }
                     }
                   }
                 }
               };

  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_results_by_individual {
  my $self = shift;
  my $individual_data = shift;
  my $ar_data = shift;
  my $results_by_individual = {};

  foreach my $individual (keys %$individual_data) {
    foreach my $gene_id (keys %{$individual_data->{$individual}}) {
      foreach my $transcript_id (keys %{$individual_data->{$individual}->{$gene_id}}) {
        foreach my $allelic_requirement (keys %{$ar_data->{$gene_id}}) {
          my $zyg2var = $individual_data->{$individual}->{$gene_id}->{$transcript_id};
          my $filtered_zyg2var = $self->zyg2var_filtered_by_allelic_requirement_rule($allelic_requirement, $zyg2var);
          if (defined $filtered_zyg2var) {
            $results_by_individual->{$individual}->{$gene_id}->{$allelic_requirement}->{$transcript_id} = $filtered_zyg2var;
          }
        }
      }
    }
  }
  return $results_by_individual;
}

=head2 get_highest_frequency

  Arg [1]    : String $frequency
  Example    : my $highest_frequency = get_highest_frequency('AA=0.0007289,AFR=0,AMR=0.0014,EA=0.005595,EAS=0,EUR=0.006,SAS=0,gnomAD=0.002923,gnomAD_AFR=0.0005175');
  Description: Extract the highest frequency from the frequency string.
  Returntype : Float $highest_frequency
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub get_highest_frequency {
  my $frequencies = shift;
  my $highest_frequency = 0;
  foreach my $frequency_annotation (split(',', $frequencies)) {
    my $frequency = (split('=', $frequency_annotation))[-1];
    if ($frequency > $highest_frequency) {
      $highest_frequency = $frequency;
    }
  }
  return $highest_frequency;
}

=head2 store_population_names

  Arg [1]    : String $frequencies
  Example    : $self->store_population_names('AA=0.0007289,AFR=0,AMR=0.0014,EA=0.005595,EAS=0,EUR=0.006,SAS=0,gnomAD=0.002923,gnomAD_AFR=0.0005175');
  Description: Extract population names from the frequency string and store population names in internal cache.
               Population names are used for creating the header in the HTML output file.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
sub store_population_names {
  my $self = shift;
  my $frequencies = shift;
  foreach my $frequency_annotation (split(',', $frequencies)) {
    my $population_name = (split('=', $frequency_annotation))[0];
    $self->{population_names}->{$population_name} = 1;
  }
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($vf) = @{parse_line({format => 'vcf', minimal => 1}, $line)};
  return $vf;
}

sub stats_html_head {
    my $html =<<SHTML;
<html>
<head>
  <title>VEP summary</title>
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
  <style>
    a.inactive {
      color: grey;
      pointer-events:none;
    }
  </style>
</head>
<body>
SHTML
  return $html;
}

sub stats_html_tail {
  my $script =<<SHTML;
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
  <script type="text/javascript" src="https://www.google.com/jsapi"></script>
  <script>
    \$( "input[type=checkbox]" ).on( "click", function(){
      if (\$('.target').is(':checked')) {
        \$( "div.not_canonical" ).hide();
        \$("a.not_canonical").addClass("inactive");
      } else {
        \$( "div.not_canonical" ).show();
        \$("a.not_canonical").removeClass("inactive");
      }
    } );
  \$(document).ready(function(){
    \$('[data-toggle="tooltip"]').tooltip(); 
  });
  </script>
SHTML
  return "\n</div>\n$script\n</body>\n</html>\n";
}

1;
