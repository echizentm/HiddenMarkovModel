package HiddenMarkovModel;
use strict;
use warnings;
use Data::Dumper;

use constant {
    DEFAULT_ITERATIONS => 10,
};

sub new {
    my ($class, %args) = @_;
    my $self = \%args;
    return bless($self, $class);
}

sub estimate {
    my ($self, %args) = @_;
    return unless ($args{data}   and ref($args{data})   eq 'ARRAY');
    return unless ($args{latent} and ref($args{latent}) eq 'ARRAY');
    return unless (@{$args{data}} == @{$args{latent}});

    $self->_initialize(
        data   => $args{data},
        latent => $args{latent},
    );

    for (1..DEFAULT_ITERATIONS) {
        $self->_expectation;
        $self->_maximization;
    }

    $self->_fill_latent;
    return [ map { $self->{latent_labels}[$_ - 1]; } @{$self->{latent}} ];
}

sub _initialize {
    my ($self, %args) = @_;

    $self->{size}   = @{$args{data}};
    $self->{data}   = [];
    $self->{latent} = [];

    my %data_dic;
    for (@{$args{data}}) {
        unless ($data_dic{$_}) {
            $data_dic{$_} = scalar(keys %data_dic) + 1;
        }
        push(@{$self->{data}}, $data_dic{$_});
    }
    my %latent_dic;
    for (@{$args{latent}}) {
        unless ($_) {
            push(@{$self->{latent}}, 0);
            next;
        }
        unless ($latent_dic{$_}) {
            $latent_dic{$_} = scalar(keys %latent_dic) + 1;
        }
        push(@{$self->{latent}}, $latent_dic{$_});
    }

    $self->{latent_labels} = [ sort { $latent_dic{$a} <=> $latent_dic{$b} } keys %latent_dic ];
    $self->{data_size}     = keys %data_dic;
    $self->{latent_size}   = keys %latent_dic;

    $self->{transition_vector} = []; # pi
    $self->{transition_matrix} = []; # A
    $self->{prob_given_latent} = []; # phi
    for (1..$self->{latent_size}) {
        push(@{$self->{transition_vector}}, 1 / $self->{latent_size});

        my @latent_vector;
        for (1..$self->{latent_size}) { push(@latent_vector, 1 / $self->{latent_size}); }
        push(@{$self->{transition_matrix}}, \@latent_vector);

        my @data_vector;
        for (1..$self->{data_size}) { push(@data_vector, 1 / $self->{data_size}); }
        push(@{$self->{prob_given_latent}}, \@data_vector);
    }
}

sub _expectation {
    my ($self, %args) = @_;

    $self->{alpha} = [];
    $self->{beta}  = [];

    for my $k (0..($self->{latent_size} - 1)) {
        $self->{alpha}[0][$k] = $self->{transition_vector}[$k]
                      * $self->{prob_given_latent}[$k][$self->{data}[0] - 1];
        $self->{beta}[$self->{size} - 1][$k] = 1;
    }

    for my $n (1..($self->{size} - 1)) {
        my $cn = $self->{size} - $n - 1;

        for my $k1 (0..($self->{latent_size} - 1)) {
            my $sum_a = 0;
            my $sum_b = 0;
            for my $k2 (0..($self->{latent_size} - 1)) {
                $sum_a += $self->{alpha}[$n - 1][$k2]
                        * $self->{transition_matrix}[$k2][$k1];
                $sum_b += $self->{beta}[$cn + 1][$k2]
                        * $self->{transition_matrix}[$k1][$k2]
                        * $self->{prob_given_latent}[$k2][$self->{data}[$cn + 1] - 1];
            }
            $self->{alpha}[$n][$k1] = $self->{prob_given_latent}[$k1][$self->{data}[$n] - 1]
                                    * $sum_a;
            $self->{beta}[$cn][$k1] = $sum_b;
        }
    }
}

sub _maximization {
    my ($self, %args) = @_;

    my $den_pi = 0;
    for my $k (0..($self->{latent_size} - 1)) {
        $den_pi += $self->_gamma(0, $k);
    }

    for my $k1 (0..($self->{latent_size} - 1)) {
        $self->{transition_vector}[$k1] = $self->_gamma(0, $k1) / $den_pi;

        for my $k2 (0..($self->{latent_size} - 1)) {
            my $num_a = 0;
            my $den_a = 0;
            for my $n (1..($self->{size} - 1)) {
                $num_a += $self->_xi($n - 1, $n, $k1, $k2);
                for my $k3 (0..($self->{latent_size} - 1)) {
                    $den_a += $self->_xi($n - 1, $n, $k1, $k3);
                }
            }
            $self->{transition_matrix}[$k1][$k2] = $num_a / $den_a;
        }

        for my $i (0..($self->{data_size} - 1)) {
            my $num_phi = 0; 
            my $den_phi = 0;
            for my $n (0..($self->{size} - 1)) {
                $num_phi += $self->_gamma($n, $k1) if ($self->{data}[$n] == ($i + 1));
                $den_phi += $self->_gamma($n, $k1);
            }
            $self->{prob_given_latent}[$k1][$i] = $num_phi / $den_phi;
        }
    }
}

sub _fill_latent {
    my ($self, %args) = @_;

    for my $n (0..($self->{size} - 1)) {
        next if ($self->{latent}[$n]);

        my $max_prob = 0;
        for my $k (0..($self->{latent_size} - 1)) {
            my $prob = $self->_gamma($n, $k)
                     * $self->{prob_given_latent}[$k][$self->{data}[$n] - 1];
            if ($prob > $max_prob) {
                $self->{latent}[$n] = $k + 1;
                $max_prob = $prob;
            }
        }
    }
}

sub _gamma {
    my ($self, $n, $k) = @_;

    if ($self->{latent}[$n]) {
        return ($self->{latent}[$n] == ($k + 1)) ? 1 : 0;
    }

    my $gamma = $self->{alpha}[$n][$k]
              * $self->{beta}[$n][$k];
    return $gamma;
}

sub _xi {
    my ($self, $n1, $n2, $k1, $k2) = @_;

    if ($self->{latent}[$n1] and $self->{latent}[$n2]) {
        return ($self->{latent}[$n1] == ($k1 + 1) and
                $self->{latent}[$n2] == ($k2 + 1)) ? 1 : 0;
    }

    my $xi = $self->{alpha}[$n1][$k1]
           * $self->{beta}[$n2][$k2]
           * $self->{transition_matrix}[$k1][$k2]
           * $self->{prob_given_latent}[$k2][$self->{data}[$n2] - 1];
    return $xi;
}

1;

__END__

