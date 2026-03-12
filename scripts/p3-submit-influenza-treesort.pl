=head1 Submit an Influenza Reassortment Analysis Job

This script submits a request to build a reassortment of influenza virus segments. The reassortment uses the phylogenetic
tree of a reference segment to infer which other segments belong with it.

=head1 Usage Synopsis

    p3-submit-influenza-treesort [options] output-path output-name

Start an influenza reassortment analysis job, producing output in the specified workspace path, using the specified name
for the base filenameof the output files.

=head2 Command-Line Options

=over 4

=item --workspace-path-prefix

Base workspace directory for relative workspace paths.

=item --workspace-upload-path

Name of workspace directory to which local files should be uploaded.

=item --overwrite

If a file to be uploaded already exists and this parameter is specified, it will be overwritten; otherwise, the script will error out.

=item --fasta

The name of a FASTA input file. The file must consist of whole segments with specifically formatted headers. This format
is described below in the L</FASTA Input Format> section.

=item --ref

The name of the reference segment. This must be one of the segments included in the FASTA file. The default is C<HA>.

=item --names

A comma-delimited list of segment names. The names must correspond to the segments included in the FASTA file. The
permissible segment names are: PB2, PB1, PA, HA, NP, NA, MP, NS. The default is C<HA,NA>. If the reference segment
name is not included, it will be added automatically.

=item method

The method to use for tree building. The options are C<local> and C<mincut>. The default is C<local>.

=item inf

Inference method to use. The options are C<FastTree> and C<IQ-Tree>. The default is C<FastTree>.

=item max-dev

The maximum deviation allowed from the standard substitution rate. The default is 2.0, which means that the maximum rate can be twice as high or
twice as low as the standard rate.

=item cutoff

The cutoff p-value for the reassortment tests. The default is 0.001 (1 percent).

=item clades-file

The name for an optional output file where clades with evidence of reassortment will be stored.

=back

The following options are used for assistance and debugging.

=over 4

=item --help

Display the command-line usage and exit.

=item --dry-run

Display the JSON submission string and exit without invoking the service or uploading files.

=back

=head2 FASTA Input Format

The reassortment service (internal name TreeSort) B<only> works with influenza nucleotide sequences.

The segment name must be surrounded by C<|> characters and can be followed by a date in the format YYYYC<->MMC<->DD.
For example:

    >A/swine/Michigan/A02635726/2021|1B.2.1|1998B|TTTPPT|HA|2021-04-23

where the segment is C<HA> and the date is C<2021-04-23> (April 23, 2021).

TreeSort can accept 3 different strain name formats within a FASTA header:

=over 4

=item *

The strain name is everything that remains after the segment name is removed.

=item *

The strain name starts with "EPI_ISL_" followed by a numeric (integer) value.

=item *

A strain name starts with A, B, C, or D followed by 3 to 5 spans of text inside C</> characters. 
For example: 

    A/swine/Iowa/A02635718/2021. 
    
Note that C<|> characters are not allowed in the strain name.

=back

=cut

use strict;
use Getopt::Long;
use Bio::KBase::AppService::Client;
use P3AuthToken;
use Data::Dumper;
use Bio::KBase::AppService::CommonSpec;
use Bio::KBase::AppService::GenomeIdSpec;
use Bio::KBase::AppService::UploadSpec;
use List::Util;

use constant SEGMENT_NAMES => { "PB2" => 1, "PB1" => 1, "PA" => 1, "HA" => 1, "NP" => 1, "NA" => 1, "MP" => 1, "NS" => 1 };


# Insure we're logged in.
my $p3token = P3AuthToken->new();
if (! $p3token->token()) {
    die "You must be logged into BV-BRC to use this script.";
}
# Get a common-specification processor, an uploader, and a reads-processor.
my $commoner = Bio::KBase::AppService::CommonSpec->new();
my $uploader = Bio::KBase::AppService::UploadSpec->new($p3token);

# Get the application service helper.
my $app_service = Bio::KBase::AppService::Client->new();

# Declare the option variables and their defaults.
my $sequences;
my $trimThreshold;
my $gapThreshold;
my $dnaFlag;
my $substitutionModel;
my $recipe = 'RAxML';
# Now we parse the options.
GetOptions($commoner->options(),
        'sequences=s@' => \$sequences,
        'trim-threshold=f' => \$trimThreshold,
        'gap-threshold=f' => \$gapThreshold,
        'dna' => \$dnaFlag,
        'substitution-model=s' => \$substitutionModel,
        'recipe=s' => \$recipe
        );
# Verify the argument count.
if (! $ARGV[0] || ! $ARGV[1]) {
    die "Too few parameters-- output path and output name are required.";
} elsif (scalar @ARGV > 2) {
    die "Too many parameters-- only output path and output name should be specified.  Found : \"" . join('", "', @ARGV) . '"';
}
# Validate the tuning parameters.
if (! RECIPES->{$recipe}) {
    die "Invalid recipe specified.";
}
if ($substitutionModel && ! SUB_MODELS->{$substitutionModel}) {
    die "Invalid substitution model.";
}
# Get the user sequence files.
my $type = ($dnaFlag ? 'feature_dna_fasta' : 'feature_protein_fasta');
my $sequenceFiles = $uploader->fix_file_list($sequences, $type);
# Add the type to each sequence file.
$sequenceFiles = [ map { { filename => $_, type => 'FASTA' } } @$sequenceFiles ];
# Compute the alphabet.
my $alphabet = ($dnaFlag ? "DNA" : "Protein");
# Handle the output path and name.
my ($outputPath, $outputFile) = $uploader->output_spec(@ARGV);
# Build the parameter structure.
my $params = {
    sequences => $sequenceFiles,
    alphabet => $alphabet,
    recipe => $recipe,
    output_path => $outputPath,
    output_file => $outputFile,
};
# Add the optionals.
if (defined $gapThreshold) {
    $params->{gap_threshold} = $gapThreshold;
}
if (defined $trimThreshold) {
    $params->{trim_threshold} = $trimThreshold;
}
if ($substitutionModel) {
    $params->{substitution_model} = $substitutionModel;
}
# Submit the job.
$commoner->submit($app_service, $uploader, $params, GeneTree => 'gene phylogeny');
