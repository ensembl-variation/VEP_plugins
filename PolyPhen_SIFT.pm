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

 PolyPhen_SIFT

=head1 SYNOPSIS

 mv PolyPhen_SIFT.pm ~/.vep/Plugins

 # Read default PolyPhen/SIFT SQLite file in $HOME/.vep
 ./vep -i variations.vcf -cache --plugin PolyPhen_SIFT

 # Read database with custom name and/or located in a custom directory
 ./vep -i variations.vcf -cache --plugin PolyPhen_SIFT,db=custom.db
 ./vep -i variations.vcf -cache --plugin PolyPhen_SIFT,dir=/some/custom/dir
 ./vep -i variations.vcf -cache --plugin PolyPhen_SIFT,db=custom.db,dir=/some/custom/dir

 # Create PolyPhen/SIFT SQLite file based on Ensembl database
 ./vep -i variations.vcf -cache --plugin PolyPhen_SIFT,create_db=1

 # Only get SIFT prediction and score
 ./vep -i variations.vcf -cache --plugin PolyPhen_SIFT,db=custom.db,polyphen=o,humdiv=o

=head1 DESCRIPTION

 A VEP plugin that retrieves PolyPhen and SIFT predictions from a
 locally constructed SQLite database. It can be used when your main
 source of VEP transcript annotation (e.g. a GFF file or GFF-based cache)
 does not contain these predictions.

 You must create a SQLite database of the predictions or point to the SQLite
 database file already created. Compatible SQLite databases based on pangenome
 data are available at http://ftp.ensembl.org/pub/current_variation/pangenomes

 You may point to the file by adding parameter `db=[file]`. If the file is not
 in `HOME/.vep`, you can also use parameter `dir=[dir]` to indicate its path.

 --plugin PolyPhen_SIFT,db=[file]
 --plugin PolyPhen_SIFT,db=[file],dir=[dir]

 To create a SQLite database using PolyPhen/SIFT data from the Ensembl database,
 you must have an active database connection (i.e. not using `--offline`) and
 add parameter `create_db=1`. This will create a SQLite file named
 `[species].PolyPhen_SIFT.db`, placed in the directory specified by the `dir`
 parameter:

 --plugin PolyPhen_SIFT,create_db=1
 --plugin PolyPhen_SIFT,create_db=1,dir=/some/specific/directory

 *** NB: this will take some hours! ***

 When creating a PolyPhen_SIFT by simply using `create_db=1`, you do not need to
 specify any parameters to load the appropriate file based on the species:

 --plugin PolyPhen_SIFT

 By default, this plugin gives SIFT score and prediction, Polyphen humvar and
 humdiv score and prediction. You can manipulate what you want using the following 
 options -

  sift      [p|s|b|o] provides SIFT prediction term, score, or both if the value is
            respectively 'p', 's', or 'b'. If the value is 'o' then do not provide SIFT
            prediction or score. Default value is 'b'.
  polyphen  [p|s|b|o] provides PolyPhen humvar prediction term, score, or both if the 
            value is respectively 'p', 's', or 'b'. If the value is 'o' then do not 
            provide PolyPhen humvar prediction or score. Default value is 'b'.
  humdiv    [p|s|b|o] provides PolyPhen humdiv prediction term, score, or both if the 
            value is respectively 'p', 's', or 'b'. If the value is 'o' then do not 
            provide PolyPhen humdiv prediction or score. Default value is 'b'.

=cut

package PolyPhen_SIFT;

use strict;
use warnings;
use DBI;
use Digest::MD5 qw(md5_hex);
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub _create_meta_table {
  my $self = shift;

  $self->{dbh}->do("CREATE TABLE meta(key, value, PRIMARY KEY (key, value))");
  my $sth = $self->{dbh}->prepare("INSERT INTO meta VALUES(?, ?)");

  my $mysql = $self->{va}->db->dbc->prepare(qq{
    SELECT meta_key, meta_value
    FROM meta m
    WHERE meta_key NOT LIKE 'schema_%'
  }, {mysql_use_result => 1});

  my ($key, $value);
  $mysql->execute();
  $mysql->bind_columns(\$key, \$value);
  $sth->execute($key, $value) while $mysql->fetch();
  $sth->finish();
  $mysql->finish();
  return 1;
}

sub _create_predictions_table {
  my $self = shift;

  $self->{dbh}->do("CREATE TABLE predictions(md5, analysis, matrix)");

  my $sth = $self->{dbh}->prepare("INSERT INTO predictions VALUES(?, ?, ?)");

  my $mysql = $self->{va}->db->dbc->prepare(qq{
    SELECT m.translation_md5, a.value, p.prediction_matrix
    FROM translation_md5 m, attrib a, protein_function_predictions p
    WHERE m.translation_md5_id = p.translation_md5_id
    AND p.analysis_attrib_id = a.attrib_id
    AND a.value IN ('sift', 'polyphen_humdiv', 'polyphen_humvar')
  }, {mysql_use_result => 1});

  my ($md5, $attrib, $matrix);
  $mysql->execute();
  $mysql->bind_columns(\$md5, \$attrib, \$matrix);
  $sth->execute($md5, $attrib, $matrix) while $mysql->fetch();
  $sth->finish();
  $mysql->finish();

  $self->{dbh}->do("CREATE INDEX md5_idx ON predictions(md5)");
  return 1;
}

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  my $species = $self->config->{species} || 'homo_sapiens';
  $self->{sift} = $param_hash->{sift} || 'b';
  $self->{polyphen} = $param_hash->{polyphen} || 'b';
  $self->{humdiv} = $param_hash->{humdiv} || 'b';

  my $dir = $param_hash->{dir} || $self->{config}->{dir};
  my $db = $param_hash->{db} || $dir.'/'.$species.'.PolyPhen_SIFT.db';

  # create DB?
  if($param_hash->{create_db}) {
    die("ERROR: DB file $db already exists - remove and re-run to overwrite\n") if -e $db;

    $self->{dbh} = DBI->connect("dbi:SQLite:dbname=$db","","");
    $self->{va} ||= Bio::EnsEMBL::Registry->get_adaptor($species, 'variation', 'variation');
    print "Creating meta table with assemblies list...\n";
    $self->_create_meta_table();

    print "Creating predictions table...\n";
    $self->_create_predictions_table();
  }

  die("ERROR: DB file $db not found - you need to download or create it first, see documentation in plugin file\n") unless -e $db;

  $self->{initial_pid} = $$;
  $self->{db_file} = $db;

  $self->{dbh} ||= DBI->connect("dbi:SQLite:dbname=$db","","");
  $self->{get_sth} = $self->{dbh}->prepare("SELECT md5, analysis, matrix FROM predictions WHERE md5 = ?");

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my ($self) = shift;
  my %header;

  $header{PolyPhen_humdiv_score}  = 'PolyPhen humdiv score from PolyPhen_SIFT plugin' if $self->{humdiv} eq 'b' || $self->{humdiv} eq 's';
  $header{PolyPhen_humdiv_pred}   = 'PolyPhen humdiv prediction from PolyPhen_SIFT plugin' if $self->{humdiv} eq 'b' || $self->{humdiv} eq 'p';
  $header{PolyPhen_humvar_score}  = 'PolyPhen humvar score from PolyPhen_SIFT plugin' if $self->{polyphen} eq 'b' || $self->{polyphen} eq 's';
  $header{PolyPhen_humvar_pred}   = 'PolyPhen humvar prediction from PolyPhen_SIFT plugin' if $self->{polyphen} eq 'b' || $self->{polyphen} eq 'p';
  $header{SIFT_score}             = 'SIFT score from PolyPhen_SIFT plugin' if $self->{sift} eq 'b' || $self->{sift} eq 's';
  $header{SIFT_pred}              = 'SIFT prediction from PolyPhen_SIFT plugin' if $self->{sift} eq 'b' || $self->{sift} eq 'p'; 
  
  return \%header;
}

sub run {
  my ($self, $tva) = @_;
  
  # only for missense variants
  return {} unless grep {$_->SO_term eq 'missense_variant'} @{$tva->get_all_OverlapConsequences};

  my $tr = $tva->transcript;
  my $tr_vep_cache = $tr->{_variation_effect_feature_cache} ||= {};

  ## if predictions are not available for both tools in the cache, look in the SQLite database
  unless(exists($tr_vep_cache->{protein_function_predictions}) &&
     $tva->sift_prediction() && $tva->polyphen_prediction()
   ){

    # get peptide
    unless($tr_vep_cache->{peptide}) {
      my $translation = $tr->translate;
      return {} unless $translation;
      $tr_vep_cache->{peptide} = $translation->seq;
    }

    # get data, indexed on md5 of peptide sequence
    my $md5 = md5_hex($tr_vep_cache->{peptide});

    my $data = $self->fetch_from_cache($md5);

    unless($data) {

      # forked, reconnect to DB
      if($$ != $self->{initial_pid}) {
        $self->{dbh} = DBI->connect("dbi:SQLite:dbname=".$self->{db_file},"","");
        $self->{get_sth} = $self->{dbh}->prepare("SELECT md5, analysis, matrix FROM predictions WHERE md5 = ?");

        # set this so only do once per fork
        $self->{initial_pid} = $$;
      }

      $self->{get_sth}->execute($md5);

      $data = {};

      while(my $arrayref = $self->{get_sth}->fetchrow_arrayref) {
        my $analysis = $arrayref->[1];
        next unless ($analysis =~ /sift|polyphen/i);
        my ($super_analysis, $sub_analysis) = split('_', $arrayref->[1]);

        $data->{$analysis} = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
          -translation_md5    => $arrayref->[0],
          -analysis           => $super_analysis,
          -sub_analysis       => $sub_analysis,
          -matrix             => $arrayref->[2]
        );
      }

      $self->add_to_cache($md5, $data);
    }

    $tr_vep_cache->{protein_function_predictions} = $data;
  }

  my $return = {};

  foreach my $tool_string(qw(SIFT PolyPhen_humdiv PolyPhen_humvar)) {
    my ($tool, $analysis) = split('_', $tool_string);
    my $lc_tool = lc($tool);

    my $check_mode = defined $analysis && $analysis eq 'humdiv' ? $analysis : $lc_tool;
    my $pred_meth  = $lc_tool.'_prediction' if ($self->{$check_mode} eq 'b' || $self->{$check_mode} eq 'p');
    my $score_meth = $lc_tool.'_score' if ($self->{$check_mode} eq 'b' || $self->{$check_mode} eq 's');

    if (defined $pred_meth) {
      my $pred = $tva->$pred_meth($analysis);

      if($pred) {
        $pred =~ s/\s+/\_/g;
        $pred =~ s/\_\-\_/\_/g;
        $return->{$tool_string.'_pred'} = $pred;
      }
    }

    if (defined $score_meth) {
      my $score = $tva->$score_meth($analysis);
      $return->{$tool_string.'_score'} = $score if defined($score);
    }
  }

  return $return;
}

sub fetch_from_cache {
  my $self = shift;
  my $md5 = shift;

  my $cache = $self->{_cache} ||= [];

  my ($data) = map {$_->{data}} grep {$_->{md5} eq $md5} @$cache;
  return $data;
}

sub add_to_cache {
  my $self = shift;
  my $md5 = shift;
  my $data = shift;

  my $cache = $self->{_cache} ||= [];
  push @$cache, {md5 => $md5, data => $data};

  shift @$cache while scalar @$cache > 50;
}

1;
