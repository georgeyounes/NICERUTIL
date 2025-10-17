"""
Build a .txt style GTI expression and .fits GTI file either
excluding TSTART-to-TSTOP times or only including them
It can be run from command line tool "buildnicergti"
"""

import pandas as pd
import os
import argparse

def build_heasoft_time_filter(df: pd.DataFrame, *, outside: bool = True) -> str:
    """
    Build a HEASoft-style time selection expression.

    Parameters
    ----------
    df : DataFrame with columns ['tstart', 'tstop']
        Each row is a closed interval [tstart, tstop].
    outside : bool
        If True (default), build expression to KEEP data OUTSIDE these intervals
        (exclude the flares). If False, build expression to KEEP data INSIDE.

    Returns
    -------
    expr : str
        A parenthesized HEASoft expression string.
    """
    if not {'tstart', 'tstop'}.issubset(df.columns):
        raise ValueError("df must have columns ['tstart','tstop'].")

    if df.empty:
        return ""

    work = df[['tstart', 'tstop']].sort_values('tstart').reset_index(drop=True)
    n = len(work)

    if outside:
        # replicate your original parenthesization pattern
        pieces = []
        for jj in range(n):
            start = work.loc[jj, 'tstart']
            stop = work.loc[jj, 'tstop']

            if jj == 0 and n == 1:
                pieces.append(f"(TIME<{start}).or.(TIME>{stop})")
            elif jj == 0:
                pieces.append(f"(TIME<{start}).or.((TIME>{stop})")
            elif jj == n - 1:
                pieces.append(f".and.(TIME<{start})).or.(TIME>{stop})")
            else:
                pieces.append(f".and.(TIME<{start})).or.((TIME>{stop})")
        expr = "".join(pieces)
    else:
        # inside intervals only
        expr = ".or.".join(
            f"((TIME>{row.tstart}).and.(TIME<{row.tstop}))" for _, row in work.iterrows()
        )

    return expr


def write_gti(df: pd.DataFrame, mkffile: str, output_prefix: str, *,
              outside: bool = True) -> str:
    """
    Build the HEASoft GTI expression from df and write to `<output_prefix>_gti.txt` and <output_prefix>_gti.fits`
    Returns the expression string
    :param df: DataFrame with columns ['tstart', 'tstop']
    :type df: pd.DataFrame
    :param mkffile: NICER mkffile (could also provide an event file since all filtering is TIME related)
    :type mkffile: str
    :param output_prefix: name of .txt and .fits GTI files ({output_prefix}_gti.txt and {output_prefix}_gti.fits)
    :type output_prefix: str
    :param outside: Keep only outside ['tstart', 'tstop'] (default) or only inside
    :type outside: bool
    :return: Expression string
    :rtype: str
    """
    expr = build_heasoft_time_filter(df, outside=outside)
    outpath = f"{output_prefix}_gti.txt"
    with open(outpath, "w") as f:
        f.write(expr)
    command = f"maketime {mkffile} {output_prefix}_gti.fits @{output_prefix}_gti.txt anything anything TIME no clobber=yes"
    os.system(command)
    return expr


def main():
    """
    CLI entry point for `buildnicergti`

    Examples:
      buildnicergti ni6533052001.mkf --output-prefix out --tstart 100 200 --tend 150 220
      buildnicergti ni6533052001.mkf --output-prefix out --tstart 100 --tend 150 --inside
    """
    parser = argparse.ArgumentParser(
        description="Build a HEASoft GTI expression (.txt) and GTI file (.fits) from TSTART/TEND intervals."
    )
    parser.add_argument("mkffile", help="NICER mkffile (or event file).")
    parser.add_argument("tstart", nargs="+", type=float, help="One or more start times (floats).")
    parser.add_argument("tend", nargs="+", type=float, help="One or more end times (floats).")
    parser.add_argument("-op", "--output-prefix", required=True, help="Prefix for output files (creates <prefix>_gti.txt and <prefix>_gti.fits).")
    parser.add_argument("-in", "--inside", action="store_true", help="Keep only INSIDE the [tstart, tstop] windows (default is to keep OUTSIDE).")
    args = parser.parse_args()

    tstarts = args.tstart
    tends = args.tend

    if len(tstarts) != len(tends):
        raise SystemExit("ERROR: --tstart and --tend must have the same number of values.")

    df = pd.DataFrame({"tstart": tstarts, "tstop": tends})

    # Optional sanity check: ensure tstart < tstop for all intervals
    if (df["tstop"] <= df["tstart"]).any():
        bad = df[df["tstop"] <= df["tstart"]]
        raise SystemExit(f"ERROR: Found intervals with tstop <= tstart:\n{bad}")

    outside = not args.inside
    expr = write_gti(df, mkffile=args.mkffile, output_prefix=args.output_prefix, outside=outside)
    print(f"Wrote {expr} to {args.output_prefix}_gti.txt file\n"
          f"Created {args.output_prefix}_gti.fits file correspoding to the requested intervals")


if __name__ == "__main__":
    main()
