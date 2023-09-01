#!/usr/bin/env perl
$host = shift;
$instance = shift;
$arg = shift;

#### random sleep, rand() can be a fraction of second
select(undef,undef,undef,rand());

if ($arg) {
  @ids = split(/,/, $arg);
}
else {
  while(1) {
    if (opendir(DDIR, "combined_relevant_proteins.faa-seq")) { 
      @ids = grep {/^\d+$/} readdir(DDIR);
      last;
    }
    else {
      sleep(1);
    }
  }
}

foreach $id (@ids) {

  next unless (-e "combined_relevant_proteins.faa-seq/$id");
  next if (-e "combined_relevant_proteins.faa-seq/$id.lock");
  $cmd = `touch combined_relevant_proteins.faa-seq/$id.lock`;

  if (50) {
    $cmd = `blastp -outfmt 6 -db ./combined_relevant_proteins.faa.1011376 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query combined_relevant_proteins.faa-seq/$id -out combined_relevant_proteins.faa-bl/$id`;
    $cmd =                         `../bin/cdhit/psi-cd-hit/psi-cd-hit.pl -J parse_blout_multi combined_relevant_proteins.faa-bl/$id -c 0.3 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 0 >> combined_relevant_proteins.faa-blm/$host.$instance`;
  }
  elsif (1) {
    $cmd = `blastp -outfmt 6 -db ./combined_relevant_proteins.faa.1011376 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query combined_relevant_proteins.faa-seq/$id | ../bin/cdhit/psi-cd-hit/psi-cd-hit.pl -J parse_blout combined_relevant_proteins.faa-bl/$id -c 0.3 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 1`;
  }
  else {
    $cmd = `blastp -outfmt 6 -db ./combined_relevant_proteins.faa.1011376 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query combined_relevant_proteins.faa-seq/$id -out combined_relevant_proteins.faa-bl/$id`;
    $cmd =                         `../bin/cdhit/psi-cd-hit/psi-cd-hit.pl -J parse_blout combined_relevant_proteins.faa-bl/$id -c 0.3 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 0`;
  }
  $cmd = `rm -f  combined_relevant_proteins.faa-seq/$id`;
  $cmd = `rm -f  combined_relevant_proteins.faa-seq/$id.lock`;
}

($tu, $ts, $cu, $cs) = times();
$tt = $tu + $ts + $cu + $cs;
$cmd = `echo $tt >> combined_relevant_proteins.faa-seq/host.$host.$instance.cpu`;

