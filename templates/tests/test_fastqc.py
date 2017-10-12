import os
import shutil
import pytest
import templates.fastqc as fq


@pytest.fixture
def temp_env():
    # Create a temporary directory and move there
    tmp_dir = ".temp"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    os.chdir(tmp_dir)
    yield
    # Exit and remove temporary directory
    os.chdir(os.pardir)
    shutil.rmtree(tmp_dir)


def test_adapters(temp_env):

    err = []

    adapters_file = "../tests/test_data/adapters.fasta"

    out = fq.convert_adatpers(adapters_file)

    with open("fastqc_adapters.tab") as fh:
        res = fh.read()

    if not isinstance(out, str):
        err.append("convert_adatpers shoud return a string: {}".format(
            str(out)))

    if res != "TruSeq_Universal_Adapter\\tAATGATACGGCGACCACCGAGATCTACA" \
              "CTCTTTCCCTACACGACGCTCTTCCGATCT\\nTruSeq_Adapter_Index 1" \
              "\\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCG" \
              "TCTTCTGCTTG\\n":
        err.append("Mismatch in output file content: {}".format(res))

    assert err == []


def test_no_adapters(temp_env):

    adapters_file = "non_existent"

    res = fq.convert_adatpers(adapters_file)

    assert res is None


def test_fastqc_run_no_adapters(temp_env):

    args = [
        ["../tests/test_data/sample_1.fastq.gz",
         "../tests/test_data/sample_2.fastq.gz"],
        None,
        "2"
    ]

    fq.main(*args)

    with open("fastq_status") as fh:
        res = fh.read()

    assert 1


def test_fastqc_run_adapters(temp_env):

    args = [
        ["../tests/test_data/sample_1.fastq.gz",
         "../tests/test_data/sample_2.fastq.gz"],
        "../tests/test_data/adapters.fasta",
        "2"
    ]

    fq.main(*args)

    with open("fastq_status") as fh:
        res = fh.read()

    assert 1
