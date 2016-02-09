#!/usr/bin/perl
use Bio::SeqIO;
use Bio::Seq;

$subtype = $ARGV[0]; # N1 or N2 or WA

if ($subtype =~ "N1")
{ $in  = Bio::SeqIO->new(-file => "A_Brisbane_59_2007_H1N1_DNA.fasta" , '-format' => 'Fasta'); }
if ($subtype =~ "N2")
{$in  = Bio::SeqIO->new(-file => "A_Brisbane_10_2007_H3N2_DNA.fasta" , '-format' => 'Fasta'); }

if ($subtype =~ "WA")
{$in  = Bio::SeqIO->new(-file => "WA.fa" , '-format' => 'Fasta'); }


while(my $seq1 = $in->next_seq) {
($segid) = $seq1->id =~m/^seg(\d+)/ ;
$segment_seq[$segid]= $seq1->seq;

if ($subtype =~ "WA" && $segid==6) { 
                $orfstart=19; 
                $segment_seq[$segid]=substr( $seq1->seq, $orfstart, 2000); 
                }

}

if ($subtype =~ "N1")
{ open F, "map-N1seg6-3BEQ" || die "no mapping file\n" ; }
if ($subtype =~ "N2")
{ open F, "map-N2seg6-1ING" || die "no mapping file\n" ; }

if ($subtype =~ "WA")
{ open F, "map-WAseg6-3BEQ" || die "no mapping file\n" ; }


while(<F>) {
chomp;
( $aa_refseq, $pos_refseq, $aa_pdb, $pos_pdb) = split(" ",$_);
$mapping{$pos_refseq} = $pos_pdb;
}
close F;


#print "$seg\n";

while(<STDIN>) {

chomp; $mutin=$_;

($number, $mutcode) = $mutin =~ m/^(.*) (s.*)$/;

($number, $seg, $wtnuc, $pos, $mutnuc) = $mutin =~ m/^(.*) s(\d+)(\w)(\d+)(\w)$/;
$pos++; # since residue id's are zero-based previously
#print "$seg $wtnuc $pos $mutnuc\n";

$pos-=$orfstart;

$old_nt = $segment_seq[$seg];
$new_nt = $segment_seq[$seg];

$oldwt = substr($new_nt, $pos-1, 1);
unless ($oldwt=~$wtnuc) { die "wt nucleotide mismatch, seq is $oldwt \n"; }
substr($new_nt, $pos-1, 1) = $mutnuc;

#print "$old_nt\n\n$new_nt\n";

$ntS_new = Bio::Seq->new(-seq=> $new_nt);
$protS_new = $ntS_new->translate;

$ntS_old = Bio::Seq->new(-seq=> $old_nt);           
$protS_old = $ntS_old->translate;

$prot_old = $protS_old->seq;
$prot_new = $protS_new->seq;

$flag=0;
for($i=0;$i<=length($prot_old);$i++) {

$aaold = substr($prot_old,$i,1) ;
$aanew = substr($prot_new,$i,1);

unless ( $aaold eq $aanew ) 
{
$p = $mapping{$i};
if (!$p) { $p="--$i--";}
print "$number $aaold$p$aanew $mutcode\n";
$flag=1;
}

}

if (!$flag)  { print "$number syn $mutcode\n"; } 

} # end while <stdin>
