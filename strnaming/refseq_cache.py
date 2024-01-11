#
# Copyright (C) 2023 Jerry Hoogenboom
#
# This file is part of STRNaming, an algorithm for generating simple,
# informative names for sequenced STR alleles in a standardised and
# automated manner.
#
# STRNaming is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STRNaming is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with STRNaming.  If not, see <http://www.gnu.org/licenses/>.
#
import gzip
import re
import time
import urllib.request
import urllib.error

from pathlib import Path

# Valid chromosome position numbers are within these ranges.
CHR_RANGES = {
    "1": (10001, 248946422),
    "2": (10001, 242183529),
    "3": (10001, 198235559),
    "4": (10001, 190204555),
    "5": (10001, 181478259),
    "6": (60001, 170745979),
    "7": (10001, 159335973),
    "8": (60001, 145078636),
    "9": (10001, 138334717),
    "10": (10001, 133787422),
    "11": (60001, 135076622),
    "12": (10001, 133265309),
    "13": (16000001, 114354328),
    "14": (16000001, 106883718),
    "15": (17000001, 101981189),
    "16": (10001, 90228345),
    "17": (60001, 83247441),
    "18": (10001, 80263285),
    "19": (60001, 58607616),
    "20": (60001, 64334167),
    "21": (5010001, 46699983),
    "22": (10510001, 50808468),
    "X": (10001, 156030895),
    "Y": (10001, 57217415),
    "M": (1, 16569)
}

# Pattern to match a chromosome number.
PAT_CHR = re.compile("^(?:[Cc][Hh][Rr])?([1-9XYM]|1[0-9]|2[0-2])$")

# Constants for cache file size and location.
CHUNK_SIZE = 12288  # Will always fit a 4 kiB block, assuming a worst-case 2/3 compression ratio.
BUILTINDIR = Path(__file__).parent / "data"
CACHEDIR = Path.home() / ".strnaming-cache"
CHUNKFILENAME = "refseq-{chromosome}-{chunk}.txt.gz"

# Constants for automatic downloading.
ENSEMBL_URL = "https://rest.ensembl.org/sequence/region/human/chr{chromosome}:{start}..{end}?content-type=text/plain"
DOWNLOAD_ATTEMPTS = 5  # Maximum number of attempts to download sequence from Ensembl.
DOWNLOAD_INTERVAL = 2  # Seconds to wait (can be float) after a failed download attempt.

# Constants for cache write access locking.
LOCK_FILE = CACHEDIR / "cache.lock"
LOCK_TIMEOUT = 60  # Maximum number of seconds to try acquiring the lock.
LOCK_INTERVAL = 0.25 # Seconds to wait (can be float) after a failed lock attempt.

# Cache write access locking mechanism.
class LockFile:
    def __init__(self):
        self.lock = None
    def __enter__(self):
        self.acquire()
    def __exit__(self, exc_type, exc_value, traceback):
        self.release()
    def __del__(self):
        self.release()
    def acquire(self):
        start_time = time.time()
        while True:
            try:
                self.lock = LOCK_FILE.open("x")
                return
            except FileExistsError:
                if time.time() - start_time > LOCK_TIMEOUT:
                    raise ValueError(
                        "Failed to acquire lock within %i seconds; please remove lock file "
                        "manually if this problem persists (%s)" % (LOCK_TIMEOUT, LOCK_FILE))
                else:
                    time.sleep(LOCK_INTERVAL)
    def release(self):
        if self.lock:
            try:
                self.lock.close()
            except:
                pass
            try:
                LOCK_FILE.unlink()
            except:
                pass
            self.lock = None


def _load_from_ensembl(chromosome, start, end):
    """Download and return the given portion of chromosome sequence."""
    for attempt in range(DOWNLOAD_ATTEMPTS):
        try:
            time.sleep(0.5)
            return urllib.request.urlopen(
                ENSEMBL_URL.format(chromosome="MT" if chromosome == "M" else chromosome,
                    start=start, end=end)).read().decode("UTF-8")
        except urllib.error.HTTPError as e:
            if e.code == 429 and "Retry-After" in e.headers:
                # We are being rate limited by the server.
                time.sleep(float(e.headers["Retry-After"]))
            else:
                time.sleep(DOWNLOAD_INTERVAL)
            last_exception = e
        except Exception as e:
            time.sleep(DOWNLOAD_INTERVAL)
            last_exception = e
    raise last_exception


def _get_chunk(chromosome, chunk, skip, length):
    """Read and return a chunked portion of reference sequence."""
    chunkfilename = CHUNKFILENAME.format(chromosome=chromosome, chunk=chunk)
    chunkfile = BUILTINDIR / chunkfilename
    try:
        with gzip.open(str(chunkfile), "rt") as f:
            f.seek(skip)
            return f.read(length)
    except FileNotFoundError:
        pass  # This chunk is not bundled with STRNaming.
    chunkfile = CACHEDIR / chunkfilename
    chunkfilepath = str(chunkfile)
    try:
        with gzip.open(chunkfilepath, "rt") as f:
            f.seek(skip)
            return f.read(length)
    except FileNotFoundError:
        # This chunk is not in the cache; try to get it.
        CACHEDIR.mkdir(parents=True, exist_ok=True)
        chunk_offset = chunk * CHUNK_SIZE + 1
        try:
            with LockFile():
                # Try reading the chunk file again, it could have been
                # created while we were waiting to acquire the lock.
                try:
                    with gzip.open(chunkfilepath, "rt") as f:
                        f.seek(skip)
                        return f.read(length)
                except FileNotFoundError:
                    pass  # Nope, still no file. We'll make it.
                seq = _load_from_ensembl(chromosome, chunk_offset,
                    min(chunk_offset + CHUNK_SIZE - 1, CHR_RANGES[chromosome][1]))
                tempfile = CACHEDIR / (chunkfilename + ".tmp")
                with gzip.open(str(tempfile), "wt") as f:
                    f.write(seq)
                try:
                    tempfile.replace(chunkfile)
                except IOError:
                    if not chunkfile.is_file():
                        raise
                return seq[skip : skip + length]
        except:
            start = chunk_offset + skip
            raise ValueError(
                "Sequence for chr{chromosome}:{start}..{end} is unavailable and could not be "
                "downloaded automatically. Please run 'strnaming "
                "refseq-cache chr{chromosome}:{start}..{end}' on a system with internet "
                "access and manually transfer the downloaded files to '{cachedir}'.".format(
                chromosome=chromosome, start=start, end=start + length - 1, cachedir=CACHEDIR))


def get_refseq(chromosome, start, end):
    """
    Get a portion of reference sequence. The end position is inclusive.
    """
    match = PAT_CHR.match(chromosome)
    if match is None:
        raise ValueError("Invalid chromosome number '%s'" % chromosome)
    chromosome = match.group(1)
    if not CHR_RANGES[chromosome][0] <= start <= end <= CHR_RANGES[chromosome][1]:
        raise ValueError(
            "Invalid refseq range: %i..%i is not within chromosome %s limits: %i..%i" % (
            start, end, chromosome, CHR_RANGES[chromosome][0], CHR_RANGES[chromosome][1]))
    start_chunk = (start - 1) // CHUNK_SIZE
    end_chunk = (end - 1) // CHUNK_SIZE
    seq = ""
    for chunk in range(start_chunk, end_chunk + 1):
        chunk_offset = chunk * CHUNK_SIZE
        skip = max(0, start - chunk_offset - 1)
        length = min(CHUNK_SIZE, end - chunk_offset) - skip
        seq += _get_chunk(chromosome, chunk, skip, length)
    return seq

