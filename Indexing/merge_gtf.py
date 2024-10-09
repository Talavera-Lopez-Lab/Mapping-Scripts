from dataclasses import dataclass, fields
from functools import reduce
import argparse

@dataclass
class GTF_File_Header:
    description: str
    provider: str
    contact: str
    format: str
    date: str

    def __add__(self, other):
        add_header_attr = lambda file1, file2, att: f"{getattr(file1, att)}; {getattr(file2, att)}"
        return GTF_File_Header(
            description = add_header_attr(self, other, "description"),
            provider = self.provider,
            contact = self.contact,
            format = self.format,
            date = add_header_attr(self, other, "date"),
        )

@dataclass
class GTF_File:
    header: GTF_File_Header
    body: list[str]

    def __add__(self, other):
        header = self.header + other.header
        body = self.body + other.body
        return GTF_File(header, body)

    def print_to_file(self, filename: str):
        if filename == None:
            filename = "concatenated.gtf"
        with open(filename, 'w') as file:
            for field in fields(self.header):
                value = getattr(self.header, field.name)
                file.write(f"##{value}\n")
            for line in self.body:
                file.write(f"{line}")

def filepath_to_gtf_file(filepath: str) -> GTF_File:
    with open(filepath, "r") as file:
        file_list = list(file)
        header = GTF_File_Header(
            description=file_list[0][2:-1],
            provider=file_list[1][2:-1],
            contact=file_list[2][2:-1],
            format=file_list[3][2:-1],
            date=file_list[4][2:-1],
        )
        body = file_list[5:]
        return GTF_File(header=header, body=body)

def main(args):
    gtf_files = list(map(filepath_to_gtf_file, args.filepaths))
    output_file = reduce(lambda a, b: a + b, gtf_files)
    output_file.print_to_file(filename=args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filepaths", nargs="+", help="Files to be concatenated")
    parser.add_argument("-o", "--output", help="Name of the output file")
    args = parser.parse_args()
    main(args)



