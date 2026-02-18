from pyliftover import LiftOver
import os

# ── Chain file paths ───────────────────────────────────────────────────────────
# hg19→hg38: downloaded on-demand by pyliftover from UCSC
# hg38→T2T:  bundled in the repo (hg38ToHs1.over.chain.gz)

_HERE = os.path.dirname(os.path.abspath(__file__))
_HG38_TO_T2T_CHAIN = os.path.join(_HERE, "hg38ToHs1.over.chain.gz")

# Cache LiftOver objects so chain files are only parsed once per session
_lo_cache = {}

def _get_lo(from_db, to_db=None, chain_file=None):
    key = chain_file or f"{from_db}->{to_db}"
    if key not in _lo_cache:
        if chain_file:
            _lo_cache[key] = LiftOver(chain_file)
        else:
            _lo_cache[key] = LiftOver(from_db, to_db)
    return _lo_cache[key]


def _liftover_bed_content(lo, bed_content):
    """
    Apply a pyliftover LiftOver object to BED text content.
    Returns (converted_str, unmapped_str).
    """
    converted_lines = []
    unmapped_lines  = []

    for line in bed_content.splitlines():
        line = line.lstrip('\ufeff')

        if not line.strip() or line.startswith(('#', 'track', 'browser')):
            converted_lines.append(line)
            continue

        parts = line.strip().split('\t')
        if len(parts) < 3:
            parts = [p.strip('"') for p in line.strip().split(',')]
        if len(parts) < 3:
            parts = line.strip().split()
        if len(parts) < 3:
            unmapped_lines.append(line + "\t# Invalid Format")
            continue

        chrom = parts[0]
        if not chrom.startswith('chr'):
            chrom = f"chr{chrom}"
        if chrom == 'chrMT':
            chrom = 'chrM'

        try:
            start = int(parts[1])
            end   = int(parts[2])
        except ValueError:
            unmapped_lines.append(line + "\t# Invalid Coordinates")
            continue

        rest = parts[3:] if len(parts) > 3 else []

        new_start_list = lo.convert_coordinate(chrom, start)
        new_end_list   = lo.convert_coordinate(chrom, end)

        if not new_start_list or not new_end_list:
            unmapped_lines.append(line + "\t# Mapping Failed")
            continue

        new_chrom = new_start_list[0][0]
        new_start = new_start_list[0][1]
        new_end   = new_end_list[0][1]

        if new_chrom != new_end_list[0][0]:
            unmapped_lines.append(line + "\t# Split Mapping")
            continue

        if new_start > new_end:
            new_start, new_end = new_end, new_start

        new_line = f"{new_chrom}\t{new_start}\t{new_end}"
        if rest:
            new_line += "\t" + "\t".join(rest)
        converted_lines.append(new_line)

    return "\n".join(converted_lines), "\n".join(unmapped_lines)


def convert_bed_hg19_to_hg38(bed_content):
    """Convert hg19 BED text to hg38. Returns (converted_str, unmapped_str)."""
    lo = _get_lo('hg19', 'hg38')
    return _liftover_bed_content(lo, bed_content)


def liftover_single_region(chrom, start, end, from_build, to_build):
    """
    Liftover a single genomic region between builds.

    Supported conversions:
        hg38 → t2t   (uses bundled hg38ToHs1.over.chain.gz)
        hg19 → hg38  (pyliftover downloads chain automatically)
        hg19 → t2t   (two-step: hg19→hg38→T2T)

    Returns (new_chrom, new_start, new_end) or None if unmapped.
    """
    from_build = from_build.lower()
    to_build   = to_build.lower()

    if from_build == to_build:
        return chrom, start, end

    # Ensure chrom has 'chr' prefix (pyliftover requires it)
    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"

    if from_build in ('hg38', 'grch38') and to_build in ('t2t', 'hs1', 'chm13'):
        lo = _get_lo(None, None, chain_file=_HG38_TO_T2T_CHAIN)
        return _convert_region(lo, chrom, start, end)

    if from_build in ('hg19', 'grch37') and to_build in ('hg38', 'grch38'):
        lo = _get_lo('hg19', 'hg38')
        return _convert_region(lo, chrom, start, end)

    if from_build in ('hg19', 'grch37') and to_build in ('t2t', 'hs1', 'chm13'):
        # Two-step via hg38
        lo1 = _get_lo('hg19', 'hg38')
        mid = _convert_region(lo1, chrom, start, end)
        if mid is None:
            return None
        lo2 = _get_lo(None, None, chain_file=_HG38_TO_T2T_CHAIN)
        return _convert_region(lo2, mid[0], mid[1], mid[2])

    raise ValueError(f"Unsupported liftover: {from_build} → {to_build}")


def _convert_region(lo, chrom, start, end):
    """Convert one region with a LiftOver object. Returns (chrom, start, end) or None."""
    sl = lo.convert_coordinate(chrom, start)
    el = lo.convert_coordinate(chrom, end)
    if not sl or not el:
        return None
    new_chrom = sl[0][0]
    new_start = sl[0][1]
    new_end   = el[0][1]
    if new_chrom != el[0][0]:
        return None
    if new_start > new_end:
        new_start, new_end = new_end, new_start
    return new_chrom, new_start, new_end


def t2t_chain_available():
    """Return True if the bundled hg38→T2T chain file exists."""
    return os.path.isfile(_HG38_TO_T2T_CHAIN)
