"""
Abstract lib with common operations for nextflow templates
"""

import bz2
import gzip
import zipfile


copen = {
    "gz": gzip.open,
    "bz2": bz2.open,
    "zip": zipfile.ZipFile
}


def guess_file_compression(file_path):
    """

    Parameters
    ----------
    file_path

    Returns
    -------

    """

    file_type = None

    magic_dict = {
        b'\x1f\x8b\x08': "gz",
        b'\x42\x5a\x68': "bz2",
        b'\x50\x4b\x03\x04': "zip"
    }

    max_len = max(len(x) for x in magic_dict)

    with open(file_path, "rb") as f:
        file_start = f.read(max_len)

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return file_type

    return file_type
