#!/usr/bin/perl -w

use Bio::Seq;

$seq_obj = Bio::Seq->new(-seq        => "aaaatgggggggggggccccgtt",
                         -display_id => "#12345",
                         -desc       => "example 1",
                         -alphabet   => "dna" );

print $seq_obj->seq();