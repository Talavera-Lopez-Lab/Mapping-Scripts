import argparse

def main(args):
    for filepath in args.filepaths:
        print(filepath)
    print(args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filepaths", nargs="+")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    main(args)
