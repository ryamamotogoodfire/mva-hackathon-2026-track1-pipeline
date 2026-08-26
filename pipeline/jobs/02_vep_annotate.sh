#!/bin/bash
# MVA Track 1 step 2: download VEP cache + annotation databases, run an in-job
# sanity slice, then annotate the full proband VCF with Ensembl VEP 112 (GRCh38).
#
# Required environment:
#   DATA_DIR  Local working directory. Holds vcf/ from step 1 and receives
#             vep_cache/, dbs/, vep_out/.
# Required tools: curl, git, tabix + bgzip (htslib), gzip, perl, and Ensembl VEP
# release 112 on PATH (or as /opt/vep/src/ensembl-vep/vep inside the public
# ensemblorg/ensembl-vep:release_112.0 container image).
set -euo pipefail

: "${DATA_DIR:?set DATA_DIR to a local working directory}"

ART="$DATA_DIR/mva-track1"
CACHE="$ART/vep_cache"
DBS="$ART/dbs"
OUTDIR="$ART/vep_out"
VCF="$ART/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz"
mkdir -p "$CACHE" "$DBS" "$OUTDIR"

# Total size of a remote file, or empty when unknown:
# 1) HEAD content-length, 2) 1-byte GET content-range.
remote_size () {
  local s
  s=$(curl -fsIL --retry 3 "$1" | awk 'BEGIN{IGNORECASE=1}/^content-length:/{v=$2} END{sub(/\r/,"",v); print v}') || true
  if [ -z "$s" ]; then
    s=$(curl -fsSL --retry 3 -r 0-0 -D - -o /dev/null "$1" | awk 'BEGIN{IGNORECASE=1}/^content-range:/{v=$3} END{sub(/.*\//,"",v); sub(/\r/,"",v); print v}') || true
  fi
  printf '%s' "$s"
}

dl () {  # dl URL dest [expected_size]
  local expect="${3:-}"
  if [ -s "$2" ]; then
    local have; have=$(stat -c%s "$2" 2>/dev/null || echo 0)
    if [ -n "$expect" ] && [ "$have" = "$expect" ]; then echo "[02] exists, skip: $2 ($have bytes)"; return 0; fi
    local want; want=$(remote_size "$1" || true)
    if [ -n "$want" ] && [ "$want" = "$have" ]; then echo "[02] exists, skip: $2 ($have bytes)"; return 0; fi
    echo "[02] partial/stale/unverifiable file ($have vs ${expect:-$want}): redownloading $2"
    rm -f "$2"
  fi
  echo "[02] downloading $1"
  curl -fSL --retry 5 --retry-all-errors -sS -o "$2.tmp" "$1" && mv "$2.tmp" "$2"
}

# Segmented parallel download for multi-GB files (single-stream Ensembl FTP is slow).
dl_fast () {  # dl_fast URL dest [nseg] [expected_size]
  local url="$1" dest="$2" nseg="${3:-16}" expect="${4:-}"
  if [ -s "$dest" ]; then
    local have; have=$(stat -c%s "$dest")
    if [ -n "$expect" ] && [ "$have" = "$expect" ]; then echo "[02] exists, skip: $dest ($have bytes)"; return 0; fi
    local size; size=$(remote_size "$url" || true)
    if [ -n "$size" ] && [ "$size" = "$have" ]; then echo "[02] exists, skip: $dest"; return 0; fi
    echo "[02] partial/stale/unverifiable file: redownloading $dest"
    rm -f "$dest"
  fi
  local size; size=$(remote_size "$url" || true)
  if [ -z "$size" ] || [ "$size" -lt 100000000 ]; then
    echo "[02] plain download (size=${size:-unknown}): $url"
    dl "$url" "$dest" "$expect"; return $?
  fi
  echo "[02] segmented download: $url ($size bytes, $nseg parts)"
  rm -f "${dest}".part*
  local chunk=$(( (size + nseg - 1) / nseg )) pids=() rc=0
  for i in $(seq 0 $((nseg-1))); do
    local s=$((i*chunk)) e=$(( (i+1)*chunk - 1 ))
    [ "$e" -ge "$size" ] && e=$((size-1))
    ( for a in 1 2 3 4 5; do
        curl -fSL --retry 8 --retry-all-errors -sS -r "${s}-${e}" -o "${dest}.part$i" "$url" && break || sleep 10
      done; echo "[02] part $i finished"; ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || rc=1; done
  [ $rc -ne 0 ] && { echo "[02] FAILED $url"; return 1; }
  : > "$dest.tmp"
  for i in $(seq 0 $((nseg-1))); do cat "${dest}.part${i}" >> "$dest.tmp"; done
  mv "$dest.tmp" "$dest"; rm -f "${dest}".part*
  local got; got=$(stat -c%s "$dest")
  echo "[02] done $dest ($got of $size bytes)"
  [ "$got" = "$size" ]
}

PLUGINS="$DBS/plugins"
# plugins live under $DBS/plugins: loftee (grch38 branch, pinned) + AlphaMissense.pm
if [ ! -f "$PLUGINS/LoF.pm" ]; then
  git clone -q https://github.com/konradjk/loftee "$PLUGINS"
  (cd "$PLUGINS" && git checkout -q a46b502a68c812c8ae0c5a5721c0603fe81cae8d && rm -rf .git)
fi
if [ ! -f "$PLUGINS/AlphaMissense.pm" ]; then
  curl -fsSL -o "$PLUGINS/AlphaMissense.pm" https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/112/AlphaMissense.pm
fi
export PERL5LIB="$PLUGINS:${PERL5LIB:-}"
export PATH="$PLUGINS:$PLUGINS/maxEntScan:${PATH}"
ls "$PLUGINS" | head -6

# The large Ensembl cache tar is only needed to create the extracted cache tree.
if [ ! -d "$CACHE/homo_sapiens/112_GRCh38" ]; then
  dl_fast "https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz" "$CACHE/homo_sapiens_vep_112_GRCh38.tar.gz" 24
  echo "[02] extracting VEP cache"
  tar -xzf "$CACHE/homo_sapiens_vep_112_GRCh38.tar.gz" -C "$CACHE"
else
  echo "[02] extracted VEP cache present, skipping tar download"
fi
dl_fast "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" "$DBS/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" 8 881964081
dl "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" "$DBS/clinvar.vcf.gz" 193420805
CV="$DBS/clinvar.vcf.gz"
dl "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi" "$CV.tbi" 610582
dl_fast "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz" "$DBS/AlphaMissense_hg38.tsv.gz" 4 642961469
if [ ! -s "$DBS/AlphaMissense_hg38.tsv.gz.tbi" ]; then
  echo "[02] building tabix index for AlphaMissense"
  tabix -f -s 1 -b 2 -e 2 "$DBS/AlphaMissense_hg38.tsv.gz"
fi
dl "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz" "$DBS/human_ancestor.fa.gz" 885231320
dl "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai" "$DBS/human_ancestor.fa.gz.fai" 736
dl "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi" "$DBS/human_ancestor.fa.gz.gzi" 764520
dl_fast "https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "$DBS/gerp_conservation_scores.homo_sapiens.GRCh38.bw" 8 12617579354
ls "$CACHE"

if [ ! -s "$DBS/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
  echo "[02] gunzipping reference fasta"
  gzip -dc "$DBS/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" > "$DBS/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fi
FASTA="$DBS/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

VEP_BIN=$(command -v vep || true)
if [ -z "$VEP_BIN" ]; then VEP_BIN=/opt/vep/src/ensembl-vep/vep; fi
"$VEP_BIN" --help 2>&1 | head -3 || true

NPROC=$(nproc)
FORKS=$(( NPROC / 2 ))
if [ "$FORKS" -gt 48 ]; then FORKS=48; fi
echo "[02] nproc=$NPROC forks=$FORKS"

VEP_COMMON=(
  --offline --cache --dir_cache "$CACHE" --cache_version 112 --assembly GRCh38
  --dir_plugins "$PLUGINS"
  --species homo_sapiens
  --fasta "$FASTA"
  --format vcf
  --tab
  --force_overwrite
  --symbol --mane --canonical --numbers --hgvs
  --sift b --polyphen b
  --af --af_1kg --af_gnomade --af_gnomadg --max_af
  --variant_class --check_existing
  --regulatory
  --plugin "LoF,loftee_path:${PLUGINS}/,human_ancestor_fa:${DBS}/human_ancestor.fa.gz,gerp_bigwig:${DBS}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
  --plugin "AlphaMissense,file=${DBS}/AlphaMissense_hg38.tsv.gz,transcript_match=1"
  --custom "file=${CV},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNDN%CLNREVSTAT"
  --fork "$FORKS"
)

echo "[02] === sanity slice ==="
SANITY_VCF="$OUTDIR/sanity.vcf"
SANITY_OUT="$OUTDIR/sanity_out.tab"
gzip -dc "$VCF" | awk 'BEGIN{n=0} /^#/{print; next} {n++; if(n<=30000) print}' > "$SANITY_VCF"
"$VEP_BIN" "${VEP_COMMON[@]}" --input_file "$SANITY_VCF" --output_file "$SANITY_OUT" 2>&1 | tail -20
ls -la "$SANITY_OUT"

perl - <<'EOF'
use strict; use warnings;
my $san = "$ENV{DATA_DIR}/mva-track1/vep_out/sanity_out.tab";
open(my $fh, '<', $san) or die "cannot open $san: $!";
my (@cols, %nonempty, %cons);
my $n = 0;
while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^##/;
    if ($line =~ /^#/) { my $h = $line; $h =~ s/^#//; @cols = split /\t/, $h; next; }
    my @row = split /\t/, $line, -1;
    $n++;
    for my $i (0..$#cols) {
        my $v = $row[$i] // '';
        $nonempty{$cols[$i]}++ if $v ne '-' && $v ne '';
    }
    my %r = map { $cols[$_] => ($row[$_] // '') } 0..$#cols;
    $cons{$r{Consequence}}++;
}
print "rows: $n\n";
print "columns: ", join(",", @cols), "\n";
my @need = qw(Uploaded_variation Location Allele Gene Feature Consequence IMPACT SYMBOL CANONICAL SIFT PolyPhen AF gnomADe_AF gnomADg_AF MAX_AF CLIN_SIG ClinVar_CLNSIG LoF LoF_filter LoF_flags am_class am_pathogenicity HGVSc HGVSp EXON INTRON);
my %have = map { $_ => 1 } @cols;
my @missing = grep { !$have{$_} } @need;
print "missing columns: @missing\n" if @missing;
for my $c (@need) { printf "  %-24s %d (%.3f)\n", $c, $nonempty{$c}//0, ($nonempty{$c}//0)/($n||1); }
my @top = sort { $cons{$b} <=> $cons{$a} } keys %cons;
print "top consequences: ", join(", ", map { "$_=$cons{$_}" } @top[0..($#top<9?$#top:9)]), "\n";
my @errors;
push @errors, "missing required columns: @missing" if @missing;
push @errors, "too few rows: $n" if $n < 25000;
for my $c ("AF","SIFT","LoF","HGVSc") { push @errors, "column never populated: $c" if $have{$c} && !($nonempty{$c}//0); }
# am_class/am_pathogenicity are positionally-validated by header only: sanity slices may contain no missense variants at all
push @errors, "AF populated on <30% of sanity variants" if @cols && (($nonempty{AF}//0)/($n||1)) < 0.30;
if (@errors) { print "SANITY FAIL: @errors\n"; exit 1; }
print "SANITY OK\n";
EOF

echo "[02] === full annotation ==="
"$VEP_BIN" "${VEP_COMMON[@]}" \
  --input_file "$VCF" \
  --output_file "$OUTDIR/proband_vep.tab.bgz" --compress_output bgzip \
  \
  2>&1 | tail -20

echo "[02] full annotation rows:"
zcat "$OUTDIR/proband_vep.tab.bgz" | grep -vc '^#' || true
echo "[02] DONE"
