from collections import namedtuple
import csv
import re
from statistics import mean

import numpy as np

IDS = ["first", "first_internal", "internal", "internal_last", "last"]
Exon = namedtuple('Exon', ['locus', 'hit_index', 'hit_id'])


def get_matrices(locus_to_exon1: dict[str: Exon],
                 locus_to_exon2: dict[str: Exon]
) -> tuple[list[list[float, float]], np.ndarray]:
    """Builds matrices comparing exons' positional usage profiles.

    Args:
        locus_to_exon1: A `dict` mapping locus to `Exon` instance.
        locus_to_exon2: A `dict` mapping locus to `Exon` instance.

    Returns:
        A tuple, where tuple[0] is an index matrix and tuple[1] an ID matrix;
        i.e.,
            - tuple[0] is a 2 x n  matrix sorted on column 1, where column 1 is
              `locus_to_exon1` exons' indices and column 2 is the same exons'
              indices in `locus_to_exon2`; n is the number of common exons in
              both dicts.
            - tuple[1] is a 5 x 5  matrix, where each column 1 represents
              `locus_to_exon1` and column 2 `locus_to_exon2`.
    """
    index_matrix = []
    id_matrix = np.zeros((5, 5))

    for locus in set(locus_to_exon1.keys()) & set(locus_to_exon2.keys()):
        index_matrix.append(
            [locus_to_exon1[locus].hit_index, locus_to_exon2[locus].hit_index])

        m = IDS.index(locus_to_exon1[locus].hit_id)
        n = IDS.index(locus_to_exon2[locus].hit_id)
        id_matrix[m][n] += 1

    return sorted(index_matrix, key=lambda v: v[0]), id_matrix


# region Parsing
def read_exon(exon_path: str) -> dict[str: Exon]:
    """Reads EXON file.

    Args:
        exon_path: The EXON file path.

    Returns:
        A dict mapping locus of `Exon` instance.
    """
    locus_to_exon = {}
    with open(exon_path, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        fields = next(reader)

        for row in reader:
            row_map = {k: v for k, v in zip(fields, row)}

            locus = row_map['exon']
            hit_id = row_map['ID']
            hit_index = float(row_map['HITindex'])

            if re.match('[A-Z]', hit_id):
                hit_id = _camel_to_snake(
                    re.match(r'^\w+(?=_)', hit_id).group())

            locus_to_exon[locus] = Exon(locus, hit_index, hit_id)

    return locus_to_exon


def rep_merge(rep1_exon: str, rep2_exon: str) -> dict[str: Exon]:
    """Merges EXON files.

    Args:
        rep1_exon: The replicate EXON file path.
        rep2_exon: The replicate EXON file path.

    Returns:
        A dict mapping locus to `Exon` instance, where the exon was identified
        in both replicates with the same positional class; indices are averaged.
    """
    locus_to_exon = {}

    locus_to_rep1 = read_exon(rep1_exon)
    locus_to_rep2 = read_exon(rep2_exon)

    for locus in set(locus_to_rep1.keys()) & set(locus_to_rep2.keys()):
        rep1 = locus_to_rep1[locus]
        rep2 = locus_to_rep2[locus]

        if rep1.hit_id == rep2.hit_id:
            locus_to_exon[locus] = Exon(
                locus,
                mean([rep1.hit_index, rep2.hit_index]),
                rep1.hit_id)

    return locus_to_exon


def _camel_to_snake(camel: str) -> str:
    """Converts camel- to snake-case.

    Args:
        camel: The camel-cased string.

    Returns:
        The `camel` string in snake-case.
    """
    return "_".join(
        re.split('(?<=[a-z])(?=[A-Z])', camel)).lower()
# endregion
