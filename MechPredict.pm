
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

MechPredict  

=head1 SYNOPSIS  

 mv MechPredict.pm ~/.vep/Plugins  
 ./vep -i input.vcf --plugin MechPredict,file=mechpredict_data.tsv  

=head1 DESCRIPTION  

This is a plugin for the Ensembl Variant Effect Predictor (VEP) that annotates missense variants with predicted 
dominant-negative (DN), gain-of-function (GOF), or loss-of-function (LOF) mechanisms derived from a
Support Vector Classification (SVC) model (Badonyi et al., 2024). These probabilities do not predict a gene is 
involved in disease, they predict the most likely molecular mechanism by which deleterious variants in a given gene 
could cause disease, if it was a dominant disease gene.

Note:
- The plugin requires MechPredict_input.tsv, a pre-processed prediction dataset in TSV format. 
- The wrangled file should contain gene-level probabilities for the three mechanism categories.  
- The plugin adds the following fields to the VEP output:  
  - `MechPredict_pDN`: Probability of a **dominant-negative (DN) mechanism**  
  - `MechPredict_pGOF`: Probability of a **gain-of-function (GOF) mechanism**  
  - `MechPredict_pLOF`: Probability of a **loss-of-function (LOF) mechanism**  
  - `MechPredict_prediction`: Statement of the most likely mechanism based on empirically-derived cutoffs from Badonyi et al., 2024. 

Usage:
1. Download the Badonyi et al., 2024 raw data, available at the links below:
  - GOF: https://osf.io/h45ns
  - DN: https://osf.io/xfy38
  - LOF https://osf.io/dj4qg

2. The plugin input data can then be prepared from the raw data using:  
   ```bash
    cut --complement -f4 pdn_svm_poly_2023-07-25.tsv | awk '{print $1 " " $2 "\t" $0}' | sort >pdn_mod.tsv &&
        cut --complement -f4 pgof_svm_poly_2023-07-25.tsv | awk '{print $1 " " $2 "\t" $0}' | sort >pgof_mod.tsv &&
        cut --complement -f4 plof_svm_poly_2023-07-28.tsv | awk '{print $1 " " $2 "\t" $0}' | sort >plof_mod.tsv &&
        join -t $'\t' -1 1 -2 1 pdn_mod.tsv pgof_mod.tsv | join -t $'\t' -1 1 -2 1 - plof_mod.tsv | cut --complement -f1,5,6,8,9 | sed '1i gene uniprot_id pDN pGOF pLOF' >MechPredict_input.tsv &&
        rm pdn_mod.tsv pgof_mod.tsv plof_mod.tsv
   ```

3. VEP can be run with the MechPredict plugin as follows:  
   ```bash
   ./vep -i variations.vcf --plugin MechPredict,file=/path/to/mechpredict_data.tsv
    ```

Citation:
Badonyi M, Marsh JA (2024) Proteome-scale prediction of molecular mechanisms underlying dominant genetic diseases. 
PLoS ONE 19(8): e0307312. https://doi.org/10.1371/journal.pone.0307312

=cut

package MechPredict;

# -- Setup ---------------------------------------------------------------------

use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Data::Dumper;

# -- Initialise plugin ---------------------------------------------------------

sub new {

    my $class = shift;

    my $self = $class->SUPER::new(@_);

    my @params = @{ $self->params };

    my %params;
    foreach my $param (@params) {
        my ( $key, $value ) =
          split( '=', $param, 2 );    # Split "file=/path/to/file"
        $params{$key} = $value if defined $key and defined $value;
    }

    my $file = $params{file}
      || die "Error: No data file supplied to MechPredict.\n";

    $self->{file} = $file;

    # Read in data file
    $self->{data} = $self->read_tsv($file);

    return $self;
}

# -- Define key subroutines ----------------------------------------------------

# Define subroutine for reading in the .tsv file
sub read_tsv {

    # Retreive the plugin object
    my ( $self, $file ) = @_;

    # Declare empty hash to store data from the .tsv
    my %data;

    # Open file for reading, passing fh as the file handle for file
    # If the file cannot be opened, then exit
    open my $fh, "<", $file or die "Could not open file '$file': $!";

    # Loop over each line of the file
    while (<$fh>) {

        # Remove trailing \n chars
        chomp;

        # There are 4 cols in the .tsv file, so assign result of split to 4 usefully named variables
        my ( $gene, $uniprot_id, $pDN, $pGOF, $pLOF ) = split( "\t", $_ );

        # Store data in single hash ref per gene
        $data{$gene} = {
            uniprot_id => $uniprot_id,
            pDN        => $pDN,
            pGOF       => $pGOF,
            pLOF       => $pLOF
        };
    }
    close $fh;

    # Return the data hash
    return \%data;
}

# Defines the feature type the plugin will run on
sub feature_types {
    return ['Transcript'];
}

# Define the VEP header annotation fields
sub get_header_info {
    return {
        MechPredict_pDN =>
'Probability that the gene is associated with a dominant-negative (DN) mechanism.',
        MechPredict_pGOF =>
'Probability that the gene is associated with a gain-of-function (GOF) mechanism.',
        MechPredict_pLOF =>
'Probability that the gene is associated with a loss-of-function (LOF) mechanism.',
        MechPredict_prediction =>
'The likely disease mechanism(s) predicted for this gene, if it is disease associated.'
    };
}

# -- Main logic ----------------------------------------------------------------

# Define hash containing probability thresholds for each mechanism as per Badonyi et al., 2024
my %thresholds = (
    pdn  => 0.61,    # Probability of dominant-negative mechanism
    pgof => 0.63,    # Probability of gain-of-function mechanism
    plof => 0.64     # Probability of loss-of-function mechanism
);

sub run {

    my ( $self, $tva ) = @_;

    # Get transcript ID
    my $transcript = $tva->transcript;

    # Get gene name
    # Usually, ->get_Gene->external_name would be used, but this doesn't work in offline mode, so pull from cached value if --offline
    # _gene_symbol is a cached value in offline mode
    my $gene_name = $transcript->{_gene_symbol};

    # Return empty hash if gene name is undefined
    return {} unless $gene_name;

    # Check if the variant has a missense consequence
    # Return empty if no missense consequence is found
    return {}
      unless grep { $_->SO_term eq 'missense_variant' }
      @{ $tva->get_all_OverlapConsequences };

    # Check whether the gene_name can be found in the MechPredict prediction data
    my $gene_data = $self->{data}{$gene_name};
    return {} unless $gene_data;

    # Pull out MechPredict prediction data for gene_name
    my ( $pdn, $pgof, $plof ) = @{$gene_data}{qw(pDN pGOF pLOF)};
    
    # Compare values to thresholds and populate prediction
    # Create prediction field
    my $prediction = "";

    # Check each value against its threshold and append to prediction
    $prediction .= "gene_predicted_as_associated_with_dominant_negative_mechanism, " if $pdn  >= $thresholds{pdn};
    $prediction .= "gene_predicted_as_associated_with_gain_of_function_mechanism, " if $pgof >= $thresholds{pgof};
    $prediction .= "gene_predicted_as_associated_with_loss_of_function_mechanism, " if $plof >= $thresholds{plof};

    # Remove trailing comma and space if
    $prediction =~ s/, $//;

    # If no predictions met the threshold, assign a default message
    $prediction = "no_conclusive_mechanism_predicted" if $prediction eq "";

    # Add the data to the VEP output
    return {
        MechPredict_pDN  => $pdn,   # Probability of dominant-negative mechanism
        MechPredict_pGOF => $pgof,  # Probability of gain-of-function mechanism
        MechPredict_pLOF => $plof,  # Probability of loss-of-function mechanism
        MechPredict_prediction => $prediction
    };
}

1;
