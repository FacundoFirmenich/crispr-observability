"""Inventory a public ZIP over verified HTTP byte ranges."""
import argparse
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile


class RemoteZip(io.RawIOBase):
    def __init__(self, url, size):
        self.url, self.size, self.offset = url, size, 0
        self.tail_offset = max(0, size - 2**20)
        self.tail = self.fetch(self.tail_offset, size)

    def fetch(self, start, end):
        if start == end:
            return b""
        req = Request(self.url, headers={"Range": f"bytes={start}-{end-1}"})
        with urlopen(req, timeout=30) as response:
            if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end-1}/{self.size}":
                raise OSError(f"Range not honored: {response.status}, {response.headers.get('Content-Range')}")
            data = response.read(end - start + 1)
        if len(data) != end - start:
            raise OSError("Incomplete range response")
        return data

    def tell(self):
        return self.offset

    def seek(self, offset, whence=0):
        self.offset = (0, self.offset, self.size)[whence] + offset
        if self.offset < 0:
            raise ValueError("Negative seek")
        return self.offset

    def read(self, size=-1):
        end = self.size if size < 0 else min(self.offset + size, self.size)
        start = min(end, self.offset)
        data = self.tail[start-self.tail_offset:end-self.tail_offset] if start >= self.tail_offset else self.fetch(start, end)
        self.offset = end
        return data

    def seekable(self):
        return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    args = parser.parse_args()
    record = json.loads(args.record.read_bytes())
    entry, = record["files"]
    url = entry["links"]["self"]
    with ZipFile(RemoteZip(url, entry["size"])) as archive:
        rows = [{"member": x.filename, "bytes": x.file_size,
                 "compressed": x.compress_size, "crc32": f"{x.CRC:08x}"}
                for x in archive.infolist()]
    out = args.record.with_suffix(".zip_index.json")
    out.write_text(json.dumps({"url": url, "archive_md5": entry["checksum"],
                               "members": rows}, indent=2) + "\n", encoding="utf-8")
    print(len(rows), "members", out.name)
    for row in rows:
        if any(word in row["member"].lower() for word in ("clone", "sanger", "g2019", ".ab1", ".txt", ".xlsx", ".csv", ".pdf")):
            print(row)


if __name__ == "__main__":
    main()
