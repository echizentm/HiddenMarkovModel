#!/usr/bin/perl
use strict;
use warnings;
use Test::More;

BEGIN { use_ok('HiddenMarkovModel') };

test_01();

sub test_01 {
    note('check estimate()');

    my $hmm = HiddenMarkovModel->new();

    is($hmm->estimate(), undef, 'estimate() with no params');

    is($hmm->estimate(
        data => ['curry', 'cake', 'ramen', 'cookie'],
    ), undef, 'estimate() without latent');
    is($hmm->estimate(
        latent => ['hot', 'sweet', 'hot', 0],
    ), undef, 'estimate() without data');
    is($hmm->estimate(
        data   => ['curry', 'cake' , 'ramen', 'cookie'],
        latent => ['hot'  , 'sweet', 'hot'],
    ), undef, 'estimate() with data and latent, but inequal sizes');

    is_deeply($hmm->estimate(
        data   => ['curry', 'cake' , 'ramen', 'cookie'],
        latent => ['hot'  , 'sweet', 'hot'  , 0       ],
    ), ['hot', 'sweet', 'hot', 'sweet'], 'estimate() with valid params');
}

done_testing();

