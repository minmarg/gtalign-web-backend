#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2025 Mindaugas Margelevicius, Vilnius University
## Transform all chains of a 3D structure (model)

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Spec;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use Getopt::Long;

##constants:
my  $MYPROGNAME = basename($0);
my  $devnull = File::Spec->devnull();
my  $TARPROG = `which tar 2>$devnull`; chomp($TARPROG);
my  $GZPROG = `which gzip 2>$devnull`; chomp($GZPROG);

my  $usage = <<EOIN;

Transform all chains of a 3D structure (model).
(C) 2025 Mindaugas Margelevicius, Vilnius University

Usage:
$0 <Options>

Options:

-i <input>       Input structure file, optionally gzipped.

-o <output>      Output transformed structure file.

-m <model_num>   Input structure's model number (optional).

-r <rotation>    Rotation matrix in row-major order with comma-separated elements.

-t <translation> Translation vector with comma-separated elements.

-h               This text.

EOIN


my  $INPUT = '';
my  $OUTPUT = '';
my  $MODEL = '';
my  $ROT = '';
my  $TRN = '';
my  $Fail = 0;

my  $result = GetOptions(
               'i=s'        => \$INPUT,
               'o=s'        => \$OUTPUT,
               'm=s'        => \$MODEL,
               'r=s'        => \$ROT,
               't=s'        => \$TRN,
               'help|h'     => sub {print $usage; exit(0);}
);


do { print $usage; exit(1); } unless $result;
do { print STDERR "ERROR: Input file missing or not found.\n"; exit(1); } unless($INPUT && -f $INPUT);
do { print STDERR "ERROR: Output file not given.\n"; exit(1); } unless($OUTPUT);
do { print STDERR "ERROR: Rotation missing.\n"; exit(1); } unless($ROT);
do { print STDERR "ERROR: Translation missing.\n"; exit(1); } unless($TRN);

my  $rexcif = qr/\.cif(?:\.gz)?$/;
my  $rexpdb = qr/\.(?:pdb|ent)(?:\.gz)?$/;

unless(($INPUT =~ /$rexcif/) || ($INPUT =~ /$rexpdb/)) {
    print STDERR "ERROR: Recognized extensions for input: .cif[.gz], .(pdb|ent)[.gz]\n";
    exit(1);
}

my  @rot = split(/,/, $ROT);
my  @trn = split(/,/, $TRN);

do { print STDERR "ERROR: 9 elements expected for Rotation.\n"; exit(1); } unless($#rot == 8);
do { print STDERR "ERROR: 3 elements expected for Translation.\n"; exit(1); } unless($#trn == 2);

foreach my $re(@rot) {
    if(!looks_like_number($re)) {
        print STDERR "ERROR: Invalid Rotation element(s).\n";
        exit(1);
    }
}

foreach my $te(@trn) {
    if(!looks_like_number($te)) {
        print STDERR "ERROR: Invalid Translation element(s).\n";
        exit(1);
    }
}

if($INPUT =~ /$rexcif/) {
    exit(1) unless TransformCif($INPUT, $OUTPUT, $MODEL, \@rot, \@trn);
}
elsif($INPUT =~ /$rexpdb/) {
    exit(1) unless TransformPDB($INPUT, $OUTPUT, $MODEL, \@rot, \@trn);
}

exit(0);

## =============================================================================
##
sub TransformCif
{
    my ($input, $output, $model, $rrot, $rtrn) = @_;
    my $gzpd = ($input =~ /\.gz$/);

    if($gzpd) {
        unless(open(F, '-|', "$GZPROG -dc '$input'")) {
            print(STDERR "ERROR: Failed to open gzipped input file: $input\n");
            return 0;
        }
    } else {
        unless(open(F, "$input")) {
            print(STDERR "ERROR: Failed to open input file: $input\n");
            return 0;
        }
    }

    unless(open(O, '>', "$output")) {
        print(STDERR "ERROR: Failed to open file for writing: $output\n");
        close(F);
        return 0;
    }

    my ($f,$fmodn,$fx,$fy,$fz) = (0,-1,-1,-1,-1);

    while(<F>) {
        $f = 0 if(/^loop_\s*$/);
        if(/^_atom_site\./) {
            $fmodn = $f if /^_atom_site\.pdbx_PDB_model_num/;
            $fx = $f if /^_atom_site\.Cartn_x\s/;
            $fy = $f if /^_atom_site\.Cartn_y\s/;
            $fz = $f if /^_atom_site\.Cartn_z\s/;
            $f++;
        }
        do {print O; next} unless(/^(?:ATOM|HETATM)\s+/);
        do {print O; next} unless((0 < $fx) && (0 < $fy) && (0 < $fz));

        my @a = split(/\s+/);

        next if((length($model) > 0) && (0 < $fmodn) && ($fmodn <= $#a) && ($a[$fmodn] ne $model));
        next if(($#a < $fx) || ($#a < $fy) || ($#a < $fz));

        my $fm = $fx; $fm = $fy if $fy < $fm; $fm = $fz if $fz < $fm; $fm--;
        my ($x, $y, $z) = ($a[$fx], $a[$fy], $a[$fz]);
        my $tx = sprintf("%8.3f", $$rrot[0] * $x + $$rrot[1] * $y + $$rrot[2] * $z + $$rtrn[0]);
        my $ty = sprintf("%8.3f", $$rrot[3] * $x + $$rrot[4] * $y + $$rrot[5] * $z + $$rtrn[1]);
        my $tz = sprintf("%8.3f", $$rrot[6] * $x + $$rrot[7] * $y + $$rrot[8] * $z + $$rtrn[2]);

        s/^((?:\S+\s+){$fm}\S+)\s+\S+\s+\S+\s+\S+\s+(.*$)/$1 $tx $ty $tz $2/;
        print O;
    }

    close(O);
    close(F);

    return 1;
}

## =============================================================================
##
sub TransformPDB
{
    my ($input, $output, $model, $rrot, $rtrn) = @_;
    my $gzpd = ($input =~ /\.gz$/);

    if($gzpd) {
        unless(open(F, '-|', "$GZPROG -dc '$input'")) {
            print(STDERR "ERROR: Failed to open gzipped input file: $input\n");
            return 0;
        }
    } else {
        unless(open(F, "$input")) {
            print(STDERR "ERROR: Failed to open input file: $input\n");
            return 0;
        }
    }

    unless(open(O, '>', "$output")) {
        print(STDERR "ERROR: Failed to open file for writing: $output\n");
        close(F);
        return 0;
    }

    my ($fmodn) = ('');

    while(<F>) {
        $fmodn = $1 if /^MODEL\s+(\S+)/;
        my $endmdl = (/^ENDMDL\s/);
        my $modelsne = ((length($model) > 0) && (length($fmodn) > 0) && ($fmodn ne $model));

        next if($modelsne && !$endmdl);
        do {$fmodn = ''; next} if($modelsne && $endmdl);

        do {print O; next} unless(/^(?:ATOM|HETATM)\s+/);

        my ($x, $y, $z) = (substr($_,30,8), substr($_,38,8), substr($_,46,8));
        $x =~ s/\s//g; $y =~ s/\s//g; $z =~ s/\s//g;
        my $tx = sprintf("%8.3f", $$rrot[0] * $x + $$rrot[1] * $y + $$rrot[2] * $z + $$rtrn[0]);
        my $ty = sprintf("%8.3f", $$rrot[3] * $x + $$rrot[4] * $y + $$rrot[5] * $z + $$rtrn[1]);
        my $tz = sprintf("%8.3f", $$rrot[6] * $x + $$rrot[7] * $y + $$rrot[8] * $z + $$rtrn[2]);

        $tx = substr($tx,0,8) if length($tx) > 8;
        $ty = substr($ty,0,8) if length($ty) > 8;
        $tz = substr($tz,0,8) if length($tz) > 8;

        substr($_,30,8) = $tx; substr($_,38,8) = $ty; substr($_,46,8) = $tz;

        print O;
    }

    close(O);
    close(F);

    return 1;
}


