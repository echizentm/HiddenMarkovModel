#!/usr/bin/perl
use strict;
use warnings;
use HiddenMarkovModel;

my $hmm = HiddenMarkovModel->new();
while (my $line = <STDIN>) {
    chomp($line);
    my ($data_str, $latent_str) = split(/ : /, $line);
    my @data   = split(/ /, $data_str);
    my @latent = split(/ /, $latent_str);

    my $result = $hmm->estimate(
        data   => \@data,
        latent => \@latent,
    );
    print join(' ', @$result)."\n";
}

