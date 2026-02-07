import argparse

from .orchestrator import Orchestrator


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", required=True, help="Path to input file.")
    args = parser.parse_args()

    return Orchestrator(input_file=args.input)


if __name__ == "__main__":
    main()
