from dataclasses import dataclass, fields
from functools import reduce
import argparse

@dataclass(frozen=True)
class GTF_File_Header:
    description: str
    provider: str
    contact: str
    format: str
    date: str

    def __add__(self, other):
        add_header_attr = lambda att: f"{getattr(self, att)}; {getattr(other, att)}"
        return GTF_File_Header(
            description = add_header_attr("description"),
            provider = self.provider,
            contact = self.contact,
            format = self.format,
            date = add_header_attr("date"),
        )

@dataclass(frozen=True)
class GTF_File:
    header: GTF_File_Header
    body: list[str]

    def __add__(self, other):
        header = self.header + other.header
        body = list(set(self.body) | set(other.body))
        return GTF_File(header, body)

    @classmethod
    def from_filepath(cls, filepath: str):
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
        return cls(header=header, body=body)

    def print_to_file(self, filename: str):
        if filename == None:
            filename = "concatenated.gtf"
        with open(filename, 'w') as file:
            for field in fields(self.header):
                value = getattr(self.header, field.name)
                file.write(f"##{value}\n")
            for line in self.body:
                file.write(f"{line}")

def main(args):
    gtf_files = list(map(GTF_File.from_filepath, args.filepaths))
    output_file = reduce(lambda a, b: a + b, gtf_files)
    output_file.print_to_file(filename=args.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepaths", nargs="+", help="Files to be concatenated")
    parser.add_argument("-o", "--output", help="Name of the output file")
    args = parser.parse_args()
    main(args)



