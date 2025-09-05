"""
A module to merge nicer event files into a master event files
"""

import re
import shlex
from pathlib import Path
from typing import List, Optional, Union, Literal
import argparse

from nicerutil.convenience import *

Which = Literal["day", "night", "all"]


def find_evt_files(input_dir: Union[str, Path], which: Which) -> List[Path]:
    """
    Return a sorted list of EVT Paths found under obsid/*/xti/event_cl, filtered by 'day'/'night'/'all'.
    """
    input_dir = Path(input_dir).resolve()
    out: List[Path] = []
    for obsid, event_dir in _iter_obsid_event_dirs(input_dir):
        if which in ("day", "all"):
            p = event_dir / f"ni{obsid}_0mpu7_cl_day_bc.evt"
            if p.exists():
                out.append(p)
        if which in ("night", "all"):
            p = event_dir / f"ni{obsid}_0mpu7_cl_night_bc.evt"
            if p.exists():
                out.append(p)
    return sorted(out)


def merge_evt_gti_fpm(
    input_dir: Union[str, Path],
    which: Which,
    output_dir: Union[str, Path] = ".",
    output_prefix: Optional[str] = None,
    list_filename: Optional[Union[str, Path]] = None,
):
    """
    Three-stage merge using files discovered directly under the obsid tree.

    Stage 1: Merge EVENTS -> tmp1
    Stage 2: Replace first entry with tmp1, merge GTI -> tmp2
    Stage 3: Replace first entry with tmp2, merge FPM_SEL -> final
    """
    input_dir = Path(input_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    files = find_evt_files(input_dir, which)
    if not files:
        raise RuntimeError(f"No matching .evt files for '{which}' under {input_dir}")

    base_list = [p.resolve() for p in files]
    if output_prefix is None:
        output_prefix = f"merged_{which}"

    tmp1 = output_dir / f"{output_prefix}_tmp1_EVENTS.evt"
    tmp2 = output_dir / f"{output_prefix}_tmp2_GTI.evt"
    final = output_dir / f"{output_prefix}_final.evt"

    if list_filename is None:
        list_path = output_dir / f".ftmerge_{which}.lst"
    else:
        list_path = Path(list_filename).resolve()

    # Step 1: EVENTS -> tmp1
    _write_ftmerge_list(base_list, "EVENTS", list_path)
    execute(f"ftmerge @{shlex.quote(str(list_path))} {shlex.quote(str(tmp1))}")

    # Step 2: first = tmp1, GTI -> tmp2
    step2 = base_list.copy()
    step2[0] = tmp1
    _write_ftmerge_list(step2, "GTI", list_path)
    execute(f"ftmerge @{shlex.quote(str(list_path))} {shlex.quote(str(tmp2))}")

    # Step 3: first = tmp2, FPM_SEL -> final
    step3 = base_list.copy()
    step3[0] = tmp2
    _write_ftmerge_list(step3, "FPM_SEL", list_path)
    execute(f"ftmerge @{shlex.quote(str(list_path))} {shlex.quote(str(final))}")

    return {
        "which": which,
        "n_files": len(base_list),
        "list_file": list_path,
        "tmp1": tmp1,
        "tmp2": tmp2,
        "final": final,
    }


def _iter_obsid_event_dirs(input_dir: Path):
    """Yield (obsid, event_dir) for subdirs matching 10 digits with xti/event_cl present."""
    obsid_re = re.compile(r"^\d{10}$")
    for subdir in sorted(input_dir.iterdir()):
        if subdir.is_dir() and obsid_re.match(subdir.name):
            event_dir = subdir / "xti" / "event_cl"
            if event_dir.is_dir():
                yield subdir.name, event_dir


def _write_ftmerge_list(evt_files: List[Path], table: str, list_path: Path) -> None:
    """Write ftmerge list: absolute path per line with [TABLE] suffix."""
    with list_path.open("w") as f:
        for p in evt_files:
            f.write(f"{p.resolve()}[{table}]\n")


def main():
    parser = argparse.ArgumentParser(
        description="Merge NICER EVT files directly from an obsid tree (EVENTS -> GTI -> FPM_SEL)."
    )
    parser.add_argument("-i", "--input", default=".",
                        help="Root directory containing obsid subdirs (default: .)")
    parser.add_argument("-w", "--which", choices=["day", "night", "all"], default="all",
                        help="Which files to merge")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Directory for merged outputs and list file (default: .)")
    parser.add_argument("-p", "--prefix", default=None,
                        help="Output prefix (default: merged_<which>)")
    parser.add_argument("-l", "--list_evtfilenames", default=None,
                        help="Optional explicit path for the ftmerge list file")
    args = parser.parse_args()

    merge_evt_gti_fpm(
        input_dir=args.input,
        which=args.which,  # type: ignore[arg-type]
        output_dir=args.outdir,
        output_prefix=args.prefix,
        list_filename=args.list_evtfilenames,
    )


if __name__ == "__main__":
    main()
