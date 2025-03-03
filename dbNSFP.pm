=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

 dbNSFP

=head1 SYNOPSIS

 mv dbNSFP.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin dbNSFP,/path/to/dbNSFP.gz,col1,col2
 ./vep -i variations.vcf --plugin dbNSFP,'consequence=ALL',/path/to/dbNSFP.gz,col1,col2
 ./vep -i variations.vcf --plugin dbNSFP,'consequence=3_prime_UTR_variant&intron_variant',/path/to/dbNSFP.gz,col1,col2

=head1 DESCRIPTION

 A VEP plugin that retrieves data for missense variants from a tabix-indexed
 dbNSFP file.
 
 Please cite the dbNSFP publications alongside the VEP if you use this resource:
 - dbNSFP      https://www.ncbi.nlm.nih.gov/pubmed/21520341
 - dbNSFP v2.0 https://www.ncbi.nlm.nih.gov/pubmed/23843252
 - dbNSFP v3.0 https://www.ncbi.nlm.nih.gov/pubmed/26555599
 - dbNSFP v4   https://www.ncbi.nlm.nih.gov/pubmed/33261662
 
 You must have the 'Bio::DB::HTS' module or the tabix utility must be installed
 in your path to use this plugin.
 
 About dbNSFP data files:
   - Downoad dbNSFP files from https://sites.google.com/site/jpopgen/dbNSFP.
   - There are two distinct branches of the files provided for academic and
     commercial usage. Please use the appropriate files for your use case.
   - The file must be processed depending on dbNSFP release version and assembly
     (see commands below). We recommend using '-T' option with the sort command
     to specify a temporary directory with sufficient space.
   - The resulting file must be indexed with tabix before use by this plugin
     (see commands below).

 For release 4.7c:
 > version=4.7c
 > wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP${version}.zip
 > unzip dbNSFP${version}.zip
 > zcat dbNSFP${version}_variant.chr1.gz | head -n1 > h

 # GRCh38/hg38 data
 > zgrep -h -v ^#chr dbNSFP${version}_variant.chr* | sort -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP${version}_grch38.gz
 > tabix -s 1 -b 2 -e 2 dbNSFP${version}_grch38.gz

 # GRCh37/hg19 data
 > zgrep -h -v ^#chr dbNSFP${version}_variant.chr* | awk '$8 != "." ' | sort -k8,8 -k9,9n - | cat h - | bgzip -c > dbNSFP${version}_grch37.gz
 > tabix -s 8 -b 9 -e 9 dbNSFP${version}_grch37.gz

 When running the plugin you must list at least one column to retrieve from the
 dbNSFP file, specified as parameters to the plugin, such as:

   --plugin dbNSFP,/path/to/dbNSFP.gz,LRT_score,GERP++_RS

 You may include all columns with 'ALL'; this fetches a large amount of data per
 variant:

   --plugin dbNSFP,/path/to/dbNSFP.gz,ALL
 
 Tabix also allows the data file to be hosted on a remote server. This plugin is
 fully compatible with such a setup - simply use the URL of the remote file:

   --plugin dbNSFP,http://my.files.com/dbNSFP.gz,col1,col2

 The plugin replaces occurrences of ';' with ',' and '|' with '&'. However, some
 data field columns, e.g. 'Interpro_domain', use the replacement characters. We
 added a file with replacement logic for customising the required replacement
 of ';' and '|' in dbNSFP data columns. In addition to the default replacements
 (';' to ',' and '|' to '&') users can add customised replacements. Users can either modify
 the file 'dbNSFP_replacement_logic' in the VEP_plugins directory or provide their own
 file as second argument when calling the plugin:

   --plugin dbNSFP,/path/to/dbNSFP.gz,/path/to/dbNSFP_replacement_logic,LRT_score,GERP++_RS
 
 Note that transcript sequences referred to in dbNSFP may be out of sync with
 those in the latest release of Ensembl; this may lead to discrepancies with
 scores retrieved from other sources.

 If the dbNSFP README file is found in the same directory as the data file,
 column descriptions will be read from this and incorporated into the VEP output
 file header.

 The plugin matches rows in the tabix-indexed dbNSFP file on:
   - genomic position
   - alt allele
   - 'aaref' - reference amino acid
   - 'aaalt' - alternative amino acid

 To match only on the genomic position and the alt allele use 'pep_match=0':

   --plugin dbNSFP,/path/to/dbNSFP.gz,pep_match=0,col1,col2

 Some fields contain multiple values, one per Ensembl transcript ID.
 By default all values are returned, separated by ';' in the default VEP output format.
 To return values only for the matched Ensembl transcript ID use 'transcript_match=1'.
 This behaviour only affects transcript-specific fields; non-transcript-specific fields
 are unaffected.

   --plugin dbNSFP,/path/to/dbNSFP.gz,transcript_match=1,col1,col2

 NB 1: Using this flag may cause no value to return if the version of the Ensembl
 transcript set differs between VEP and dbNSFP.

 NB 2: MutationTaster entries are keyed on a different set of transcript IDs. Using
 the 'transcript_match' flag with any MutationTaster field selected will have no effect
 i.e. all entries are returned. Information on corresponding transcript(s) for
 MutationTaster fields can be found using http://www.mutationtaster.org/ChrPos.html

=cut

package dbNSFP;

use strict;
use warnings;
use File::Basename qw(basename);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my %INCLUDE_SO = map {$_ => 1} qw(missense_variant stop_lost stop_gained start_lost);
my %ALLOWED_PARAMS = map {$_ => 1} qw(consequence pep_match transcript_match);

# this region chosen as it should pull out a row in all assemblies
my $EXAMPLE_REGION = "1:1008170-1082927";

# these fields are ";"-separated but NOT transcript-specific OR do not correspond to Ensembl_transcriptid
# MutationTaster: Information on corresponding transcript(s) can be found by querying http://www.mutationtaster.org/ChrPos.html
my %NON_TRANSCRIPT_SPECIFIC_FIELDS = map {$_ => 1} qw(
  MutationTaster_score
  MutationTaster_converted_rankscore
  MutationTaster_pred
  MutationTaster_model
  Interpro_domain
  MutPred_Top5features
  Transcript_id_VEST3
  Transcript_var_VEST3
  VEST3_score
);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  # get dbNSFP file
  my @params = grep !/=/, @{$self->params};
  my $file = shift @params;
  my $basename = basename($file, ".gz");

  my $version;
  if ($file =~ /2\.9/) {
    $version = '2.9';
  } elsif ($file =~ /4\.0b1/) {
    # special case: different name for location column
    $version = '4.0.1';
  } elsif ($file =~ /([3-9]\.\d+(?:\.\d+)?[ac])/) {
    # prioritise versions with explicit 'a' or 'c' suffix
    $version = $1;
  } elsif ($file =~ /([3-9]\.\d+(?:\.\d+)?)/) {
    # fallback to general version, warn user
    $version = $1;
    warn "WARNING: Specific suffix 'a' or 'c' not found for version $version. Proceeding with general version.\n";
  } else {
    die "ERROR: Could not retrieve dbNSFP version from filename $file\n";
  }
  $self->{dbNSFP_version} = $version;
  $self->{basename} = $basename;

  $self->add_file($file);
  $self->get_dbNSFP_file_header($file, $EXAMPLE_REGION);

  # check if 2nd argument is a file that specifies replacement logic
  # read replacement logic 
  my $replacement_file = $params[0];
  if (defined $replacement_file && -e $replacement_file) {
    $self->add_replacement_logic($replacement_file);
    shift @params;
  } else {
    $self->add_replacement_logic();
  }

  $self->get_named_params();

  # get user-specified dbNSFP columns
  my @invalid;
  for my $col (@params) {
    if($col eq 'ALL') {
      $self->{cols} = {map {$_ => 1} @{$self->{headers}}};
      last;
    }
    
    # store columns as invalid if not found in file header
    my $valid = (grep {$_ eq $col} @{$self->{headers}});
    $self->{cols}->{$col} = $valid;
    push(@invalid, $col) unless $valid;
  }
  
  die "ERROR: no input columns were found in file header. Available columns are:\n" .
    join(",", @{$self->{headers}})."\n"
    unless defined($self->{cols}) && scalar keys %{$self->{cols}};
  
  warn "WARNING: the following columns were not found in file header: ",
       join(",", @invalid), "\n" if (@invalid);
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  if(!exists($self->{_header_info})) {

    # look for readme
    my $file_dir = $self->files->[0];

    my %rm_descs;

    # won't work for remote
    if($file_dir !~ /tp\:\/\//) {

      # get just dir
      $file_dir =~ s/\/[^\/]+$/\//;

      if(opendir DIR, $file_dir) {
        my ($readme_file) = grep {/dbnsfp.*readme\.txt/i} readdir DIR;
        closedir DIR;

        if(open RM, $file_dir.$readme_file) {
          my ($col, $reading);

          # parse dbNSFP readme
          # relevant lines look like:
          #
          # 1   column1_name: description blah blah
          #     blah blah blah
          # 2   column2_name: description blah blah
          #     blah blah blah

          while(<RM>) {
            chomp;
            s/\r$//g;

            if(/^\d+\s/) {
              $reading = 1;

              m/^\d+\s+(.+?)\:\s+(.+)/;
              $col = $1;

              $rm_descs{$col} = "(from $self->{basename}) ".$2 if $col && $2;
            }
            elsif($reading && /\w/) {
              s/^\s+//;
              $rm_descs{$col} .= ' '.$_;
            }
            else {
              $reading = 0;
            }
          }

          close RM;

          # remove multiple spaces
          $rm_descs{$_} =~ s/\s+/ /g for keys %rm_descs;
        }
      }
    }

    $self->{_header_info} = {map {$_ => $rm_descs{$_} || ($_.' from dbNSFP file')} keys %{$self->{cols}}};
  }
  
  return $self->{_header_info};
}

sub run {
  my ($self, $tva) = @_;
  
  # only for missense variants
  if ($self->{consequence} eq 'filter') {
    return {} unless grep {$INCLUDE_SO{$_->SO_term}} @{$tva->get_all_OverlapConsequences};
  }
  
  my $vf = $tva->variation_feature;
  my $tv = $tva->transcript_variation;
  
  return {} unless $vf->{start} == $vf->{end};
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT]$/;
  
  # get transcript stable ID
  my $tr_id = $tva->transcript->stable_id;

  my $data;
  my $pos;
  my $allele_string;

  my $assembly = $self->{config}->{assembly};
  my $chr = ($vf->{chr} =~ /MT/i) ? 'M' : $vf->{chr};
  foreach my $tmp_data(@{$self->get_data($chr, $vf->{start} - 1, $vf->{end})}) {
    # compare allele and transcript
    if ($assembly eq 'GRCh37') {
      if (exists $tmp_data->{'pos(1-coor)'} && $self->{dbNSFP_version} eq '2.9') {
        # for dbNSFP version 2.9.1
        $pos = $tmp_data->{'pos(1-coor)'}
      } elsif (exists $tmp_data->{'hg19_pos(1-based)'}) {
        # for dbNSFP version 3.5c indexed for hg19/(=GRCh37)
        $pos =  $tmp_data->{'hg19_pos(1-based)'}
      } else {
        die "dbNSFP file does not contain required columns (pos(1-coor) for version 2.9.1 or hg19_pos(1-based) for the other versions) to use with GRCh37\n";
      }
    } else {
      if (exists $tmp_data->{'pos(1-based)'}) {
        $pos = $tmp_data->{'pos(1-based)'}
      } elsif (exists $tmp_data->{'pos(1-coor)'} && $self->{dbNSFP_version} eq '4.0.1' ) {
        $pos = $tmp_data->{'pos(1-coor)'};
      } else {
        die "dbNSFP file does not contain required column pos(1-based) to use with GRCh38 or pos(1-coor) for dbNSFP version ".$self->{dbNSFP_version} ."\n" ;
      }
    }
    next unless
      $pos == $vf->{start} &&
      defined($tmp_data->{alt}) &&
      $tmp_data->{alt} eq $allele;

    if ($self->{pep_match}) {
      $allele_string = join('/', $tmp_data->{aaref}, $tmp_data->{aaalt});
      $allele_string =~ s/X/*/g;
      next if (!$tva->pep_allele_string() || $tva->pep_allele_string() ne $allele_string);
    }

    # make a clean copy as we're going to edit it
    %$data = %$tmp_data;

    # check and parse data if we're using transcript_match parameter
    if($self->{transcript_match}) {

      # find the "index" of this transcript
      my @tr_ids = split(';', $data->{Ensembl_transcriptid});
      my $tr_index;

      for my $i(0..$#tr_ids) {
        if($tr_ids[$i] eq $tr_id) {
          $tr_index = $i;
          last;
        }
      }

      # if transcriptid doesn't match we have to nerf transcript-specific fields
      if(!defined($tr_index)) {
        for my $key(@{$self->{transcript_specific_fields}}) {
          delete $data->{$key} if defined($data->{$key});
        }
      }

      # refine values of transcript-specific fields
      # we need to get the fields here as we're modifying $data in the loop
      my @refine_fields = grep {defined($data->{$_}) && defined($self->{cols}->{$_})} @{$self->{transcript_specific_fields}};

      foreach my $key(@refine_fields) {
	next if $data->{$key} eq '.';

        my @split = split(';', $data->{$key});

	if($tr_index > $#split) {
          warn("ERROR: Transcript index out of range for field $key\n");
	  next;
        }

	if(scalar @split != scalar @tr_ids) {
	  warn("ERROR: Number of transcript IDs does not match number of data entries for field $key\n");
	  next;
        }
	
	$data->{$key} = $split[$tr_index];
      }
    }
    last;
  }
  
  return {} unless scalar keys %$data;
  
  # get required data
  my @from = @{$self->{replacement}->{default}->{from}};
  my @to = @{$self->{replacement}->{default}->{to}};

  my %return;
  foreach my $colname (keys %{$self->{cols}}) {
    if (!$self->{cols}->{$colname}) {
      $return{$colname} = "invalid_field";
      next;
    }
    next if(!defined($data->{$colname}));
    next if($data->{$colname} eq '.');

    my @from = @{$self->{replacement}->{default}->{from}};
    my @to   = @{$self->{replacement}->{default}->{to}};
    @from    = @{$self->{replacement}->{$colname}->{from}} if (defined $self->{replacement}->{$colname});
    @to      = @{$self->{replacement}->{$colname}->{to}} if (defined $self->{replacement}->{$colname});
    for my $i (0 .. $#from) {
      $data->{$colname} =~ s/\Q$from[$i]\E/$to[$i]/g;
    }
    $return{$colname} = $data->{$colname};
  }
  
  return \%return;
}

sub parse_data {
  my ($self, $line) = @_;

  $line =~ s/\r$//g;

  my @split = split /\t/, $line;
  
  # parse data into hash of col names and values
  my %data = map {$self->{headers}->[$_] => $split[$_]} (0..(scalar @{$self->{headers}} - 1));

  return \%data;
}

sub get_start {  
  return $_[1]->{'pos(1-based)'};
}

sub get_end {
  return $_[1]->{'pos(1-based)'};
}

sub get_dbNSFP_file_header {
  my $self = shift;
  my $file = shift;
  my $region = shift;

  open HEAD, "tabix -fh $file $region 2>&1 | ";
  while(<HEAD>) {
    chomp;

    # parse header line to get field names
    if(/^\#/) {
      $_ =~ s/^\#//;
      $self->{headers} = [split];
    }

    # parse data line to identify transcript-specific fields
    else {
      next unless /\;/;
      die "ERROR: No headers found before data\n" unless defined($self->{headers});
      my $row_data = $self->parse_data($_);
      my @transcript_specific_fields;
      for my $key(keys %$row_data) {
        push @transcript_specific_fields, $key if $row_data->{$key} =~ /\;/ && !$NON_TRANSCRIPT_SPECIFIC_FIELDS{$key};
      }
      $self->{transcript_specific_fields} = \@transcript_specific_fields;
      last;
    }
  }
  close HEAD;
}

sub add_replacement_logic {
  my $self = shift;
  my $file = shift;
  $file ||= 'dbNSFP_replacement_logic';
  if (! -e $file) {
    $self->{replacement}->{default}->{from} = [';', '|'];
    $self->{replacement}->{default}->{to} = [',', '&'];
  } else {
    open FILE, $file;
    while(<FILE>) {
      chomp;
      next if /^colname/;
      my ($colname, $from, $to) = split/\s+/;
      die ("ERROR: replacement logic file requires 3 values separated by whitespace in this format: colname from to.\n")
        if(!($colname && $from && $to));
      push @{$self->{replacement}->{$colname}->{from}}, $from;
      push @{$self->{replacement}->{$colname}->{to}}, $to;
    }
    close FILE;
    die("ERROR: No default replacement logic has been specified.\n")
      if (!defined $self->{replacement}->{default});
  }

  die "ERROR: Could not read headers from $file\n"
    unless defined($self->{headers}) && scalar @{$self->{headers}};

  # check alt and Ensembl_transcriptid headers
  foreach my $h(qw(alt Ensembl_transcriptid)) {
    die "ERROR: Could not find required column $h in $file\n"
      unless grep {$_ eq $h} @{$self->{headers}};
  }
}

sub get_named_params {
  my $self = shift;

  $self->{consequence} = 'filter';
  $self->{pep_match} = 1;
  $self->{transcript_match} = 0;

  my @named_params = grep /=/, @{$self->params};
  for my $param (@named_params) {
    next if $param !~ /=/;
    my ($name, $value) = split '=', $param;

    if ($name eq "consequence") {
      # parse consequences
      if (uc $value eq 'ALL') {
        $self->{consequence} = 'ALL';
      } else {
        %INCLUDE_SO = map {$_ => 1} split/&/, $value;
      }
    } elsif ($ALLOWED_PARAMS{$name}) {
      $self->{$name} = $value;
    } else {
      die "ERROR: Invalid parameter $name\n";
    }
  }

  if ($self->{pep_match}) {
    # Check the columns for the aa are there
    foreach my $h (qw(aaalt aaref)) {
      die("ERROR: Could not find the required column $h for pep_match option in dbNSFP file\n")
        unless grep{$_ eq $h} @{$self->{headers}};
    }
  }

  if ($self->{transcript_match} && !defined($self->{transcript_specific_fields})) {
    die("ERROR: transcript_match parameter specified but transcript-specific field detection failed\n");
  }
}

1;
