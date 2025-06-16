import os, pathlib, requests, urllib.parse as _u

def download_zenodo_record(record_url_or_id, dest_dir=".",
                           overwrite=False, show_progress=True,
                           chunk_size=2**13):
    """
    Grab every file attached to a public Zenodo record (old or new API).

    Parameters
    ----------
    record_url_or_id : str | int
        Either ``https://zenodo.org/records/<recid>`` or just ``<recid>``.
    dest_dir : str          Where to put the files (directory is created).
    overwrite : bool        Re-download if an identically-sized copy exists.
    show_progress : bool    Display tqdm progress bars if the package is there.
    chunk_size : int        Streaming chunk size in bytes (default: 8192).

    Returns
    -------
    list[pathlib.Path]      Absolute paths of all downloaded files.
    """
    # ------------------------------------------------------------------ helpers
    def _record_id(x):
        """Extract the numeric recid from a URL or int-ish."""
        try:
            return int(x)
        except ValueError:
            return int([p for p in _u.urlparse(x).path.split("/") if p][-1])

    def _direct_file_url(file_entry, rec_id):
        """Return a usable download URL for one file_entry."""
        links = file_entry.get("links", {})
        return (links.get("download")              # legacy (<2023-10)
             or links.get("content")               # new API (>=2023-10)
             or links.get("self")                  # still works to stream
             or f"https://zenodo.org/records/{rec_id}/files/"
                f"{_u.quote(file_entry['key'])}?download=1")  # brute-force

    # ---------------------------------------------------------------- meta call
    recid = _record_id(record_url_or_id)
    meta  = requests.get(f"https://zenodo.org/api/records/{recid}", timeout=30)
    meta.raise_for_status()
    meta  = meta.json()
    files = meta.get("files") or meta.get("files", {}).get("entries", [])
    if not files:
        raise RuntimeError(f"Record {recid} contains no public files.")

    # ------------------------------------------------------------- dl directory
    dest = pathlib.Path(dest_dir).expanduser().resolve()
    dest.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------------- tqdm setup
    tqdm = None
    if show_progress:
        try:
            from tqdm import tqdm as _tqdm
            tqdm = _tqdm
        except ImportError:
            pass

    downloaded = []
    for f in files:
        fname   = f["key"]
        size    = f.get("size")
        url     = _direct_file_url(f, recid)
        out     = dest / fname

        if (out.exists() and not overwrite
                and (size is None or out.stat().st_size == size)):
            downloaded.append(out)
            continue

        # stream download ------------------------------------------------------
        r = requests.get(url, stream=True)
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0)) or size or None
        bar   = tqdm(total=total, unit="B", unit_scale=True, desc=fname) if tqdm else None
        with out.open("wb") as fh:
            for chunk in r.iter_content(chunk_size):
                if chunk:
                    fh.write(chunk)
                    if bar:
                        bar.update(len(chunk))
        if bar:
            bar.close()
        downloaded.append(out)

    return downloaded

# quick test ---------------------------------------------------------------
if __name__ == "__main__":
    paths = download_zenodo_record("https://zenodo.org/records/1468560",
                                   dest_dir="./zenodo_1468560")
    print("Downloaded:\n  " + "\n  ".join(map(str, paths)))
