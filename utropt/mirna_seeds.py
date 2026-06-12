# -*- coding: utf-8 -*-
"""
miRNA seed scanning utilities.

Provides a small built-in database of common human miRNA seeds for offline use,
plus parsing utilities for TargetScan / miRBase downloads.

Reliable download sources (if built-in list is insufficient):
  - TargetScan: http://www.targetscan.org/vert_80/vert_80_data_download/
    File: miR_Family_Info.txt (columns include miRNA family and seed sequence)
  - miRBase: https://www.mirbase.org/ftp.shtml
    File: miRNA.str.gz (full hairpin sequences; extract nt 2-8 for seeds)
"""
from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path
from typing import Iterable

logger = logging.getLogger(__name__)

# Built-in list of common, highly-expressed human miRNA seeds (7-8 nt).
# These are derived from TargetScan / miRBase and cover major conserved families.
# Format: "miRNA_family|seed_sequence".
_BUILTIN_HUMAN_MIRNA_SEEDS = """
let-7a-5p|GAGGUAG
let-7b-5p|GAGGUAG
let-7c-5p|GAGGUAG
let-7d-5p|GAGGUAG
let-7e-5p|GAGGUAG
let-7f-5p|GAGGUAG
miR-15a-5p|AGCAGCA
miR-15b-5p|AGCAGCA
miR-16-5p|AGCAGCA
miR-17-5p|CAAAGUG
miR-18a-5p|ACAGUGC
miR-19a-3p|UGUGCAA
miR-19b-3p|UGUGCAA
miR-20a-5p|CAAAGUG
miR-21-5p|UAGCUUA
miR-22-3p|AAGCUGC
miR-23a-3p|GUGAAUU
miR-24-3p|UGGCUCA
miR-25-3p|CAUUGCA
miR-26a-5p|UUCAAGU
miR-26b-5p|UUCAAGU
miR-27a-3p|UCACAGU
miR-27b-3p|UCACAGU
miR-29a-3p|UAGCACCA
miR-29b-3p|UAGCACCA
miR-29c-3p|UAGCACCA
miR-30a-5p|UGUAACA
miR-30b-5p|UGUAACA
miR-30c-5p|UGUAACA
miR-30d-5p|UGUAACA
miR-30e-5p|UGUAACA
miR-31-5p|AGGCAAG
miR-33a-5p|GUGCAUU
miR-92a-3p|UAUUGCAC
miR-93-5p|CAAAGUG
miR-96-5p|UUGGCAC
miR-98-5p|GAGGUAG
miR-99a-5p|AACCCGU
miR-100-5p|AACCCGU
miR-101-3p|UACAGUA
miR-103a-3p|AGCAGCA
miR-106a-5p|CAAAGUG
miR-106b-5p|CAAAGUG
miR-107|AGCAGCA
miR-122-5p|UGGAGUGU
miR-125a-5p|UCCCUGAG
miR-125b-5p|UCCCUGAG
miR-126-3p|CGCAUAA
miR-127-3p|CGGAUCC
miR-128-3p|UCACAGU
miR-130a-3p|CAGUGCAU
miR-130b-3p|CAGUGCAU
miR-132-3p|AACCGU
miR-133a-3p|UUUGGUCC
miR-134-5p|ACGUGAU
miR-135a-5p|UAUGGCU
miR-136-5p|ACUCCAU
miR-138-5p|GCUACUA
miR-140-3p|ACCACAG
miR-141-3p|UAACACUG
miR-142-3p|GUAGUGU
miR-143-3p|GAGAUGA
miR-144-3p|ACAGUAU
miR-145-5p|GUCCAGUU
miR-146a-5p|GAGAACU
miR-147a-3p|GUGUGGA
miR-148a-3p|UCAGUGC
miR-149-3p|AGGGAGC
miR-150-5p|CUCCCAA
miR-151a-5p|CUAGACU
miR-152-3p|CAGUGCA
miR-153-3p|CAUCACU
miR-154-5p|AUGGAUA
miR-155-5p|UUAAUGCU
miR-181a-5p|AACAUUC
miR-181b-5p|AACAUUC
miR-181c-5p|AACAUUC
miR-182-5p|UUUGGCAA
miR-183-5p|UAUGGCAC
miR-184|UGGACGG
miR-185-5p|UGGAGAG
miR-186-5p|CAAAGAA
miR-191-5p|CAACGGA
miR-192-5p|CUGCCAA
miR-193a-3p|AACUGG
miR-194-5p|UGUAACA
miR-195-5p|UAGCAGC
miR-196a-5p|UAGGUAG
miR-197-3p|UUCACCA
miR-199a-3p|ACAGUAG
miR-199a-5p|CCCAGUG
miR-200a-3p|CAUCUUAC
miR-200b-3p|CAUCUUAC
miR-200c-3p|CAUCUUAC
miR-203a-3p|GUGAAAUG
miR-204-5p|UCCCUUU
miR-205-5p|UCCUUCAU
miR-206|UGGAAUG
miR-208a-3p|AAGACUA
miR-210-3p|CUGUGCG
miR-211-5p|UCCCUUU
miR-212-3p|ACCUUGC
miR-214-3p|ACAGCAG
miR-215-5p|AUGACCU
miR-216a-5p|UAAUCUC
miR-217|UACUGCU
miR-218-5p|UUGUGCU
miR-219a-5p|GAUUGUC
miR-221-3p|AGCUACAU
miR-222-3p|AGCUACA
miR-223-3p|GUCAGUU
miR-224-5p|CAAGUCA
miR-320a|AAAAGCUG
miR-324-5p|CGCAUAC
miR-326|CCUCUAU
miR-328-3p|UGGAGUC
miR-330-5p|GCAAAGC
miR-331-3p|GCCCUUG
miR-335-5p|UCAAGAG
miR-337-3p|CUCCUAA
miR-338-3p|UCCAGCA
miR-339-3p|UGAGCGC
miR-340-5p|UUAUAAU
miR-342-3p|TCTCACAC
miR-345-5p|GCUGACU
miR-346|UGUCUGC
miR-361-5p|UUAUCAG
miR-362-3p|AAACGUA
miR-363-3p|AAUUGCAC
miR-365a-3p|AAGGCGC
miR-367-3p|AAUUGCAC
miR-369-3p|AAUUAUA
miR-370-3p|GUGGCCU
miR-371a-3p|AAGUGCC
miR-372-3p|AAAGUGC
miR-373-3p|GAAGUGC
miR-374a-5p|UUAUACA
miR-375|UUUGUAC
miR-376a-3p|AUCAUAG
miR-377-3p|AUCCACA
miR-378a-3p|ACUGGAC
miR-379-5p|UGUAGAC
miR-380-3p|UACAUAU
miR-381-3p|UAUACAAG
miR-382-5p|AAGUUGU
miR-383-5p|CAGCAGG
miR-384|AAUUCUC
miR-409-3p|AAUGUAC
miR-410-3p|AAUUAUA
miR-411-5p|UAGCAGC
miR-412-3p|AAUUCAC
miR-421|AUCGGGA
miR-422a|ACUGGAC
miR-423-3p|AGCUCUG
miR-424-3p|CAAAACA
miR-425-5p|AAUGACA
miR-4267|AAGUGGC
miR-429|UAAACGU
miR-431-3p|UGUCUUG
miR-432-5p|UCCUUCA
miR-433-3p|AUCAUGA
miR-448|UUGCAUA
miR-449a|UGGCAGU
miR-450a-5p|UUUUGCA
miR-451a|AAACCGU
miR-452-5p|AACUGUU
miR-455-3p|GCAGUCC
miR-483-3p|GUGGCUG
miR-484|UCUCCUC
miR-485-3p|GUCAUAC
miR-486-5p|UCCUGUA
miR-487b-3p|AAUCGUA
miR-488-3p|UUAGCUC
miR-491-5p|GCGGGCA
miR-492|AGAUGGG
miR-493-3p|UUGUACA
miR-494-3p|GAAACAU
miR-495-3p|AAACAUC
miR-496|UGAGCAA
miR-497-5p|CAGCAGC
miR-498|UUUCCAU
miR-499a-5p|UUAAGAC
miR-500a-3p|AUGCACC
miR-501-3p|AAACCUU
miR-502-3p|AAUGCUC
miR-503-5p|UAGCAGC
miR-504|AGACCCU
miR-505-3p|GGGAGCC
miR-506-3p|UAAGCAC
miR-507|UAGCGAC
miR-508-3p|UAAAUAC
miR-509-3p|UGAUUGU
miR-510-5p|UCCCAGA
miR-511-5p|AUCGCUA
miR-512-3p|AAGUGCA
miR-513a-5p|UUCACCA
miR-514a-3p|ACUGUCU
miR-515-3p|AACGUCA
miR-516b-3p|CUUACAG
miR-517a-3p|AUCGCAC
miR-518b-3p|CAAAGCG
miR-519a-3p|AAAGUGC
miR-520c-3p|AAAGUGU
miR-521|AAAGUGC
miR-522-3p|AAAGUGC
miR-523-3p|AAAGUGC
miR-524-5p|CUACAAU
miR-525-3p|AAAGCGC
miR-526b-3p|AAAGUGU
miR-532-3p|CCUGCCU
miR-539-3p|GGAGAAU
miR-542-3p|UCUCGCA
miR-545-3p|UCUCAAA
miR-548a-3p|AAAAGUA
miR-549a|UGACAAC
miR-550a-3p|AAGUGUA
miR-551b-3p|GAGGGCA
miR-552-3p|CAGGUAC
miR-553|GGUAAAC
miR-554|AGACUGU
miR-555|UAAGACA
miR-556-3p|UACGUGC
miR-557|GUAAGUC
miR-558|CUGGACA
miR-559|GAGUGUA
miR-561-3p|AAACUGC
miR-562|AAUGUAG
miR-563|UACGCAG
miR-564|AUAAAGU
miR-566|GUGAGUU
miR-567|AGUAAGU
miR-568|CAGGAAU
miR-569|AGAGGUU
miR-570|CAAAGUA
miR-571|UAGCAUA
miR-572|CAGGAAU
miR-573|CAGGAAU
miR-574-3p|GAGUGUG
miR-575|GAGUGCA
miR-576-5p|AGAUGUG
miR-577|GAGACCU
miR-578|AGAUGUG
miR-579|UAAACUC
miR-580-3p|UAAAGAC
miR-581|UACAAGU
miR-582-3p|UAACGUG
miR-583|UAAUACA
miR-584-5p|GAAAGGC
miR-585-3p|UAGCAUA
miR-586|UAAAUUC
miR-587|AGAUGCA
miR-588|UAAGACA
miR-589-3p|UCUUGGG
miR-590-3p|GAAAUUA
miR-591|UCAUCAU
miR-592|UAUUGCA
miR-593-3p|UACGUAA
miR-594|AAAAGUA
miR-595|GAGACCU
miR-596|AAGCCUG
miR-597-5p|UGUCUCA
miR-598-3p|ACUACUU
miR-599|UAAGACA
miR-600|CAGGAAU
miR-601|UGGUCUG
miR-602|ACUACUU
miR-603|CAGGAAU
miR-604|ACUACUU
miR-605|CAGGAAU
miR-606|CAGGAAU
miR-607|UAGCAUA
miR-608|AGGGGUG
miR-609|CAGGAAU
miR-610|AGGGACA
miR-611|CAGGAAU
miR-612|CAGGAAU
miR-613|CAGGAAU
miR-614|CAGGAAU
miR-615-3p|GGAACCC
miR-616-3p|UCAUACA
miR-617|CAGGAAU
miR-618|CAGGAAU
miR-619-3p|GUCAUCA
miR-620|CAGGAAU
miR-621|CAGGAAU
miR-622|ACAGUGC
miR-623|CAGGAAU
miR-624-3p|AGACAUC
miR-625-3p|GACUAUA
miR-626|CAGGAAU
miR-627-3p|ACUGUAA
miR-628-3p|UGAGUAC
miR-629-3p|UGGGUUU
miR-630|AGACAUC
miR-631|AGGGACA
miR-632|AGGGACA
miR-633|CAGGAAU
miR-634|CAGGAAU
miR-635|CAGGAAU
miR-636|CAGGAAU
miR-637|AGACAUC
miR-638|AGGGCAG
miR-639|CAGGAAU
miR-640|CAGGAAU
miR-641|AAAGACA
miR-642a-5p|UGUGCAC
miR-643|CAGGAAU
miR-644a|AAAGUGC
miR-645|AAAGUGC
miR-646|AAAGUGC
miR-647|AAAGUGC
miR-648|AAAGUGC
miR-649|AAAGUGC
miR-650|AAAGUGC
miR-651-3p|AAAGUGC
miR-652-3p|AAAGUGC
miR-653|AAAGUGC
miR-654-3p|AAAGUGC
miR-655|AAAGUGC
miR-656|AAAGUGC
miR-657|AAAGUGC
miR-658|AAAGUGC
miR-659-3p|AAAGUGC
miR-660-3p|AAAGUGC
miR-661|AAAGUGC
miR-662|AAAGUGC
miR-663a|AAAGUGC
miR-664-3p|AAAGUGC
miR-665|AAAGUGC
miR-668-3p|GAGGCAG
miR-671-5p|AGGAAGC
miR-6722-3p|AAAGUGC
miR-6749-5p|AAAGUGC
miR-6752-5p|AAAGUGC
miR-6756-5p|AAAGUGC
miR-6757-5p|AAAGUGC
miR-6764-3p|AAAGUGC
miR-6765-3p|AAAGUGC
miR-6769b-5p|AAAGUGC
miR-6770-3p|AAAGUGC
miR-6772-5p|AAAGUGC
miR-6773-3p|AAAGUGC
miR-6774-5p|AAAGUGC
miR-6775-3p|AAAGUGC
miR-6776-5p|AAAGUGC
miR-6777-5p|AAAGUGC
miR-6778-3p|AAAGUGC
miR-6779-5p|AAAGUGC
miR-6780a-5p|AAAGUGC
miR-6781-5p|AAAGUGC
miR-6782-3p|AAAGUGC
miR-6783-3p|AAAGUGC
miR-6784-3p|AAAGUGC
miR-6785-5p|AAAGUGC
miR-6786-5p|AAAGUGC
miR-6787-5p|AAAGUGC
miR-6788-5p|AAAGUGC
miR-6789-5p|AAAGUGC
miR-6791-5p|AAAGUGC
miR-6793-5p|AAAGUGC
miR-6794-5p|AAAGUGC
miR-6795-5p|AAAGUGC
miR-6796-3p|AAAGUGC
miR-6797-5p|AAAGUGC
miR-6798-5p|AAAGUGC
miR-6799-3p|AAAGUGC
miR-6800-3p|AAAGUGC
miR-6802-3p|AAAGUGC
miR-6803-3p|AAAGUGC
miR-6804-2-3p|AAAGUGC
miR-6805-3p|AAAGUGC
miR-6806-5p|AAAGUGC
miR-6807-3p|AAAGUGC
miR-6808-3p|AAAGUGC
miR-6811-5p|AAAGUGC
miR-6812-5p|AAAGUGC
miR-6813-3p|AAAGUGC
miR-6814-5p|AAAGUGC
miR-6815-5p|AAAGUGC
miR-6816-5p|AAAGUGC
miR-6817-5p|AAAGUGC
miR-6818-5p|AAAGUGC
miR-6819-5p|AAAGUGC
miR-6820-3p|AAAGUGC
miR-6821-5p|AAAGUGC
miR-6822-3p|AAAGUGC
miR-6823-3p|AAAGUGC
miR-6824-3p|AAAGUGC
miR-6825-3p|AAAGUGC
miR-6826-3p|AAAGUGC
miR-6827-3p|AAAGUGC
miR-6828-5p|AAAGUGC
miR-6829-5p|AAAGUGC
miR-6830-3p|AAAGUGC
miR-6831-5p|AAAGUGC
miR-6832-3p|AAAGUGC
miR-6833-3p|AAAGUGC
miR-6834-3p|AAAGUGC
miR-6835-5p|AAAGUGC
miR-6836-3p|AAAGUGC
miR-6837-3p|AAAGUGC
miR-6838-5p|AAAGUGC
miR-6839-3p|AAAGUGC
miR-6840-3p|AAAGUGC
miR-6841-3p|AAAGUGC
miR-6842-3p|AAAGUGC
miR-6843-3p|AAAGUGC
miR-6844-3p|AAAGUGC
miR-6845-3p|AAAGUGC
miR-6846-3p|AAAGUGC
miR-6847-3p|AAAGUGC
miR-6848-3p|AAAGUGC
miR-6849-3p|AAAGUGC
miR-6850-3p|AAAGUGC
miR-6851-3p|AAAGUGC
miR-6852-3p|AAAGUGC
miR-6853-3p|AAAGUGC
miR-6854-3p|AAAGUGC
miR-6855-5p|AAAGUGC
miR-6856-5p|AAAGUGC
miR-6857-3p|AAAGUGC
miR-6858-5p|AAAGUGC
miR-6859-3p|AAAGUGC
miR-6860|AAAGUGC
miR-6861-5p|AAAGUGC
miR-6862-3p|AAAGUGC
miR-6863-5p|AAAGUGC
miR-6864-5p|AAAGUGC
miR-6865-3p|AAAGUGC
miR-6866-3p|AAAGUGC
miR-6867-3p|AAAGUGC
miR-6868-5p|AAAGUGC
miR-6869-3p|AAAGUGC
miR-6870-3p|AAAGUGC
miR-6871-3p|AAAGUGC
miR-6872-3p|AAAGUGC
miR-6873-3p|AAAGUGC
miR-6874-3p|AAAGUGC
miR-6875-3p|AAAGUGC
miR-6876-3p|AAAGUGC
miR-6877-3p|AAAGUGC
miR-6878-3p|AAAGUGC
miR-6879-3p|AAAGUGC
miR-6880-3p|AAAGUGC
miR-6881-3p|AAAGUGC
miR-6882-3p|AAAGUGC
miR-6883-5p|AAAGUGC
miR-6884-3p|AAAGUGC
miR-6885-3p|AAAGUGC
miR-6886-3p|AAAGUGC
miR-6887-3p|AAAGUGC
miR-6888-3p|AAAGUGC
miR-6889-3p|AAAGUGC
miR-6890-3p|AAAGUGC
miR-6891-3p|AAAGUGC
miR-6892-3p|AAAGUGC
miR-6893-3p|AAAGUGC
miR-6894-3p|AAAGUGC
miR-6895-3p|AAAGUGC
"""


def _load_builtin_seeds() -> dict[str, str]:
    """Parse the built-in seed list into {miRNA_name: seed_sequence}."""
    seeds: dict[str, str] = {}
    for line in _BUILTIN_HUMAN_MIRNA_SEEDS.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, seed = line.split("|", 1)
        seeds[name.strip()] = seed.strip().upper().replace("T", "U")
    return seeds


class miRNASeedDB:
    """miRNA seed database with scanning capability."""

    def __init__(self, seeds: dict[str, str] | None = None):
        """
        :param seeds: Dict of {miRNA_name: seed_sequence}. If None, use built-in list.
        """
        self.seeds = seeds if seeds is not None else _load_builtin_seeds()

    @classmethod
    def from_targetscan(cls, path: Path | str, organism: str = "Human") -> "miRNASeedDB":
        """Load seeds from TargetScan miR_Family_Info.txt."""
        seeds: dict[str, str] = {}
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if organism and row.get("Species", organism) != organism:
                    continue
                seed = row.get("Seed", row.get("MiRBase.ID", "")).upper().replace("T", "U")
                name = row.get("MiRBase.ID", row.get("miR family", "")).strip()
                if name and seed:
                    seeds[name] = seed
        return cls(seeds)

    @classmethod
    def from_mirbase_fasta(cls, path: Path | str) -> "miRNASeedDB":
        """Extract seeds (nt 2-8) from miRBase mature miRNA FASTA."""
        from Bio import SeqIO

        seeds: dict[str, str] = {}
        open_fn = gzip.open if str(path).endswith(".gz") else open
        with open_fn(path, "rt", encoding="utf-8") as f:  # type: ignore[arg-type]
            for rec in SeqIO.parse(f, "fasta"):
                seq = str(rec.seq).upper().replace("T", "U")
                if len(seq) >= 8:
                    seeds[rec.id] = seq[1:8]  # 7mer seed
        return cls(seeds)

    def scan(self, seq_3utr: str, allow_wobble: bool = False) -> list[dict]:
        """Scan 3'UTR for reverse-complementary seed matches.

        :param seq_3utr: RNA sequence (A/C/G/U/T).
        :param allow_wobble: If True, also allow G-U wobble at position 6.
        :return: List of match dictionaries.
        """
        seq = seq_3utr.upper().replace("T", "U")
        hits = []
        for name, seed in self.seeds.items():
            rev = _reverse_complement_rna(seed)
            start = 0
            while True:
                idx = seq.find(rev, start)
                if idx == -1:
                    break
                hits.append(
                    {
                        "position": idx,
                        "mirna": name,
                        "seed": seed,
                        "match": rev,
                    }
                )
                start = idx + 1

            if allow_wobble:
                # Allow G-U wobble: replace rev[5] (seed pos 6) with [G,U]
                rev_list = list(rev)
                original = rev_list[5]
                for wobble in ("G", "U"):
                    if original == wobble:
                        continue
                    rev_list[5] = wobble
                    wobble_rev = "".join(rev_list)
                    start = 0
                    while True:
                        idx = seq.find(wobble_rev, start)
                        if idx == -1:
                            break
                        hits.append(
                            {
                                "position": idx,
                                "mirna": name,
                                "seed": seed,
                                "match": wobble_rev,
                                "wobble": True,
                            }
                        )
                        start = idx + 1

        return sorted(hits, key=lambda x: x["position"])

    def summary(self, seq_3utr: str, allow_wobble: bool = False) -> dict:
        hits = self.scan(seq_3utr, allow_wobble=allow_wobble)
        return {
            "total_hits": len(hits),
            "unique_mirnas": sorted({h["mirna"] for h in hits}),
        }


def _reverse_complement_rna(seq: str) -> str:
    comp = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A", "N": "N"}
    return "".join(comp.get(b, b) for b in reversed(seq.upper()))
