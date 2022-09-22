import os


def get_test_output_dir():
    f = "/tmp/tf-tests"
    os.makedirs(f, exist_ok=True)
    return f
