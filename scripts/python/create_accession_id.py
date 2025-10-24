"""
Module for generating stable, repeatable accession IDs using UUID version 5.
Can be used as a library or run as a standalone script for demonstration.
"""

import uuid
import random


def _get_accession_type_abbr(accession_type: str) -> str:
    """
    Returns the abbreviation for a given accession type.
    Raises ValueError if the accession type is invalid.
    """
    valid_types = {'guide': "GU", 'element': "EL", 'outer_primer': "OP", 'inner_primer': "IP"}
    if accession_type.lower() not in valid_types:
        raise ValueError(f"Invalid accession_type '{accession_type}'. Must be one of {valid_types}.")
    return valid_types[accession_type.lower()]


def _get_custom_namespace(accession_type: str) -> uuid.UUID:
    """
    Returns a custom namespace UUID based on the accession type.
    """
    NAMESPACE_DNS = uuid.NAMESPACE_DNS
    namespace = uuid.uuid5(NAMESPACE_DNS, accession_type.lower())
    return namespace


def _get_uuid(namespace: uuid.UUID, sequence: str) -> uuid.UUID:
    """
    Returns a UUID version 5 for the given namespace and sequence.
    """
    return uuid.uuid5(namespace, sequence)


def _get_short_uuid(full_uuid: uuid.UUID, sequence: str, length: int = 12) -> str:
    """
    Returns a short, reproducible substring of the UUID for use in accession IDs.
    The substring is selected using a random index seeded by the sequence.
    """
    seed = sum(ord(char) for char in sequence)
    rng = random.Random(seed)
    uuid_str = str(full_uuid).replace("-", "")
    max_start = len(uuid_str) - length
    start = rng.randint(0, max_start)
    short_uuid = uuid_str[start:start + length].upper()
    return short_uuid


def get_accession(sequence: str, accession_type: str) -> tuple[uuid.UUID, str, str]:
    """
    Generates a stable, repeatable accession ID for a given sequence and accession type.
    Args:
        sequence (str): The input sequence for which to generate the UUID.
        accession_type (str): One of 'guide', 'element', 'outer_primer', 'inner_primer'.
    Returns:
        tuple: (uuid.UUID, str, str) - UUID, accession ID, and the original sequence.
    """
    #accession_domain = "HGRM"
    accession_type_abbr = _get_accession_type_abbr(accession_type)
    accession_namespace = _get_custom_namespace(accession_type)
    sequence_uuid = _get_uuid(accession_namespace, sequence)
    accession_uuid = _get_short_uuid(sequence_uuid, sequence)
    return sequence_uuid, f"{accession_type_abbr}{accession_uuid}", sequence


def main() -> None:
    """
    Demonstrates usage of get_hgrm_accession with example sequences and types.
    """
    sequence1 = "AAAA"
    sequence2 = "AAAA"
    sequence3 = "TTTT"

    uuid1 = get_accession(sequence1, "guide")
    uuid2 = get_accession(sequence2, "guide")
    uuid3 = get_accession(sequence3, "guide")
    uuid4 = get_accession(sequence3, "element")

    print(f"UUID for 'AAAA': {uuid1}")
    print(f"UUID for 'AAAA': {uuid2}")
    print(f"UUID for 'TTTT': {uuid3}")
    print(f"UUID for 'TTTT': {uuid4}")

    # The first two UUIDs will be identical
    print(f"\nAre UUIDs for 'AAAA' identical? {uuid1 == uuid2}")
    # The third UUID will be different
    print(f"Is UUID for 'TTTT' different? {uuid1 != uuid3}")
    # The fourth UUID will be different from the third
    # because of the different accession type
    print(f"Is UUID for 'TTTT' (element) different from 'TTTT' (guide)? {uuid3 != uuid4}")


if __name__ == "__main__":
    main()
