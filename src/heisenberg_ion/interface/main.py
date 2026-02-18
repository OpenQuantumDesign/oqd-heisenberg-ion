import argparse

from .orchestrator import Orchestrator


def extract_override_arguments(args):

    override_dict = {}
    if args.override:
        for arg_kvp in args.override:
            override_dict[arg_kvp[0]] = arg_kvp[1]

    return override_dict


def main(argv=None):
    """
    main entry point of the package when called from the command line

    Args:
        argv (list[str], optional): list of command line arguments. Defaults to None.

    Returns:
        (int): exit code, equals 0 if program executes successfully
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", required=True, help="Path to input file.")
    parser.add_argument(
        "--override",
        "-o",
        required=False,
        action="append",
        nargs=2,
        metavar=("KEY", "VALUE"),
        help="Override key word arguments",
    )
    args = parser.parse_args(argv)

    override_dict = extract_override_arguments(args)

    Orchestrator(input_file=args.input, **override_dict)

    return 0


if __name__ == "__main__":
    main()
