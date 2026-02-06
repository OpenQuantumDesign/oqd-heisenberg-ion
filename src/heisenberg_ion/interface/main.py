import argparse

from .workflow import Workflow


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", required=True, help="Path to input file.")
    args = parser.parse_args()

    return Workflow(InputFile=args.input)


if __name__ == "__main__":
    main()
