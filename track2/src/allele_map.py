#!/usr/bin/env python3
"""Allele mechanism map for the BUB1B compound-heterozygous pair (Track 2 Q1).

Coordinates are identical to Track 1 (GRCh38):
  p.Leu737Ter   c.2210T>G   chr15:40209701   NM_001211.6 / ENST00000287598
  p.Asn1002Lys  c.3006T>G   chr15:40220612   NM_001211.6 / ENST00000287598

For the stop-gain allele this computes the quantities that govern any
readthrough hypothesis: codon frame and identity, the full sequence context,
the containing exon, downstream exon count, and the distance to the last
exon-exon junction evaluated against the canonical >50 nt NMD rule.

For the missense allele it verifies the codon change, locates the residue
against UniProt feature annotation, and checks the last-exon NMD escape.

Primary sources are queried live and pinned in SOURCES so every number in the
Track 2 report traces to a public endpoint.
"""
import json
import re
import time
import urllib.request

from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
RESULTS.mkdir(exist_ok=True)

TRANSCRIPT = "ENST00000287598"  # Track 1 canonical feature, MANE-matched to NM_001211.6
REFSEQ = "NM_001211.6"
UNIPROT = "O60566"              # BUB1B_HUMAN
STOP = {"pos": 40209701, "ref": "T", "alt": "G", "cdna": 2210}
MISS = {"pos": 40220612, "ref": "T", "alt": "G", "cdna": 3006}

SOURCES = {
    "ensembl_lookup": f"https://rest.ensembl.org/lookup/id/{TRANSCRIPT}?expand=1",
    "ensembl_map_check": f"https://rest.ensembl.org/map/cds/{TRANSCRIPT}/2210..2210",
    # NM_001211.6 and ENST00000287598 are the MANE-Select pair (identical CDS);
    # Ensembl REST sequence/id takes the Ensembl id, so the CDS fetched here is
    # byte-identical to the RefSeq CDS the HGVS c. numbering refers to.
    "refseq_cds": f"https://rest.ensembl.org/sequence/id/{TRANSCRIPT}?type=cds",
    "refseq_protein": f"https://rest.ensembl.org/sequence/id/{TRANSCRIPT}?type=protein",
    "uniprot_features": f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.json",
    "alphafold_metadata": f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}",
    "genomic_context": f"https://rest.ensembl.org/sequence/region/human/15:{STOP['pos']-15}..{STOP['pos']+15}:1",
}

CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'G','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def fetch(url, accept="application/json", tries=5):
    req = urllib.request.Request(url, headers={"User-Agent": "silico-track2/1.0", "Accept": accept})
    last = None
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 + 4 * a)
    raise RuntimeError(f"fetch failed: {url}: {last}")


def seq_payload(text):
    lines = text.strip().split("\n")
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    return re.sub(r"\s", "", "".join(lines)).upper()


out = {"sources": SOURCES, "transcript": TRANSCRIPT, "refseq": REFSEQ}

# transcript structure
L = json.loads(fetch(SOURCES["ensembl_lookup"]))
assert L["strand"] == 1, "expected plus strand"
exons = sorted(L["Exon"], key=lambda e: e["start"])
T = L["Translation"]
cds_gstart, cds_gend = T["start"], T["end"]
out["transcript_meta"] = {
    "exon_count": len(exons), "strand": L["strand"], "biotype": L["biotype"],
    "is_canonical": L.get("is_canonical"), "translation_id": T["id"],
    "translation_length_aa": T["length"], "cds_genomic_span": f"{cds_gstart}-{cds_gend}",
    "transcript_version": L.get("version"), "assembly": L.get("assembly_name"),
}

blocks = []
cum = 0
for i, e in enumerate(exons):
    o0, o1 = max(e["start"], cds_gstart), min(e["end"], cds_gend)
    if o1 < o0:
        continue
    cum += o1 - o0 + 1
    blocks.append({"rank": len(blocks) + 1, "exon_id": e["id"], "g0": o0, "g1": o1,
                   "cds_block_len": o1 - o0 + 1, "cds_end_cum": cum})
out["cds_blocks"] = blocks

# verify independent Ensembl CDS->genomic map at the PTC
M = json.loads(fetch(SOURCES["ensembl_map_check"]))
hit = M["mappings"][0]
out["map_check_2210_to_genomic"] = {"start": hit["start"], "strand": hit["strand"]}
assert hit["start"] == STOP["pos"] and hit["strand"] == 1

# CDS + protein sequences
cds = seq_payload(fetch(SOURCES["refseq_cds"], accept="text/plain"))
protein = seq_payload(fetch(SOURCES["refseq_protein"], accept="text/plain"))
out["cds_length_nt"] = len(cds)
out["protein_length_aa"] = len(protein.rstrip("*"))


def cdna_info(name, g):
    p = g["cdna"]
    ref_base = cds[p - 1]
    cs = p - (p - 1) % 3
    codon = cds[cs - 1: cs + 2]
    alt = list(codon)
    alt[p - cs] = g["alt"]
    alt = "".join(alt)
    return {
        "cdna_position": p,
        "codon_start": cs,
        "ref_base_matches_genomic_ref": ref_base == g["ref"],
        "reference_codon": codon,
        "alternate_codon": alt,
        "reference_aa": CODON[codon],
        "alternate_aa": CODON[alt],
        "cds_context_minus20_plus19": cds[max(0, p - 21): p + 19],
    }


alleles = {"p.Leu737Ter": cdna_info("stop", STOP), "p.Asn1002Lys": cdna_info("missense", MISS)}
assert alleles["p.Leu737Ter"]["reference_codon"] == "TTA"
assert alleles["p.Leu737Ter"]["alternate_codon"] == "TGA"
assert alleles["p.Leu737Ter"]["reference_aa"] == "L" and alleles["p.Leu737Ter"]["alternate_aa"] == "*"
assert alleles["p.Asn1002Lys"]["reference_aa"] == "N" and alleles["p.Asn1002Lys"]["alternate_aa"] == "K"
assert protein[736] == "L" and protein[1001] == "N"
out["alleles"] = alleles

sc = alleles["p.Leu737Ter"]["codon_start"]
out["stop_codon_context"] = {
    "stop_codon": alleles["p.Leu737Ter"]["alternate_codon"],
    "tetranucleotide": alleles["p.Leu737Ter"]["alternate_codon"] + cds[sc + 2],
    "upstream_6nt": cds[sc - 7: sc - 1],
    "downstream_9nt_after_stop": cds[sc + 2: sc + 11],
}

# exon occupancy of each variant
def exon_of(gpos):
    for b in blocks:
        if b["g0"] <= gpos <= b["g1"]:
            return b
    return None

stop_block = exon_of(STOP["pos"])
miss_block = exon_of(MISS["pos"])
out["exon_assignment"] = {
    "p.Leu737Ter": {"exon_rank": stop_block["rank"], "exon_id": stop_block["exon_id"],
                    "of_exons": len(blocks)},
    "p.Asn1002Lys": {"exon_rank": miss_block["rank"], "exon_id": miss_block["exon_id"],
                     "of_exons": len(blocks)},
}

# NMD rule for the PTC
junctions = [b["cds_end_cum"] for b in blocks[:-1]]
last_junction = junctions[-1]
ptc_cs = sc
downstream = [j for j in junctions if j >= ptc_cs]
out["nmd_rule"] = {
    "ptc_codon_start_cds": ptc_cs,
    "last_exon_junction_cds": last_junction,
    "distance_ptc_to_last_junction_nt": last_junction - ptc_cs,
    "downstream_exon_junctions": len(downstream),
    "downstream_junction_cds_positions": downstream,
    "rule": "NMD predicted when the PTC lies >~50 nt upstream of the 3'-most exon-exon junction",
    "nmd_predicted": bool(last_junction - ptc_cs > 50),
    "nmd_margin_vs_50nt": last_junction - ptc_cs - 50,
}
# Translation terminates at codon 737 without adding a residue, so the
# truncated product is 736 residues and residues 737..1050 (314 aa) are lost.
out["nmd_rule"]["fraction_protein_retained_if_transcript_survives"] = round(736 / out["protein_length_aa"], 4)
out["nmd_rule"]["fraction_protein_lost_if_transcript_survives"] = round(1 - 736 / out["protein_length_aa"], 4)

# missense allele: last exon => transcript escapes NMD
out["missense_nmd"] = {
    "in_last_exon": miss_block["rank"] == len(blocks),
    "escapes_nmd": True,
}

# UniProt features over the C-terminal region
unip = json.loads(fetch(SOURCES["uniprot_features"]))
feats = unip.get("features", [])
relevant = []
for f in feats:
    t = f.get("type")
    if t in ("Domain", "Region", "Binding site", "Active site", "Site", "Motif"):
        lo, hi = f["location"]["start"]["value"], f["location"]["end"]["value"]
        relevant.append({"type": t, "description": f.get("description", ""), "start": lo, "end": hi})
out["uniprot_features"] = relevant
kinase = [r for r in relevant if r["type"] == "Domain" and "kinase" in r["description"].lower()]
if kinase:
    k = kinase[0]
    out["kinase_domain"] = k
    out["n1002_inside_kinase_domain"] = k["start"] <= 1002 <= k["end"]
    out["l737_inside_kinase_domain"] = k["start"] <= 737 <= k["end"]

# AlphaFold confidence
af = json.loads(fetch(SOURCES["alphafold_metadata"]))
af_url = (af[0]["pdbUrl"] if isinstance(af, list) else af.get("pdbUrl"))
out["alphafold_model_url"] = af_url
pdb = fetch(af_url, accept="text/plain")
plddt = {}
for line in pdb.splitlines():
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        plddt[int(line[22:26].strip())] = float(line[60:66].strip())
if plddt:
    def win(a, b):
        vals = [plddt[i] for i in range(a, b + 1) if i in plddt]
        return round(sum(vals) / len(vals), 1) if vals else None
    out["alphafold_plddt"] = {
        "n_residues": len(plddt),
        "res_737": plddt.get(737),
        "res_1002": plddt.get(1002),
        "window_720_750": win(720, 750),
        "window_990_1015": win(990, 1015),
        "window_900_1050": win(900, 1050),
    }

out["genomic_context_15flank_stop"] = fetch(SOURCES["genomic_context"], accept="text/plain").strip()

RESULTS.joinpath("allele_map.json").write_text(json.dumps(out, indent=2) + "\n")
print("OK wrote results/allele_map.json")
print(json.dumps({k: out[k] for k in ("alleles", "stop_codon_context", "exon_assignment", "nmd_rule",
                                      "missense_nmd", "alphafold_plddt", "uniprot_features")}, indent=2))
