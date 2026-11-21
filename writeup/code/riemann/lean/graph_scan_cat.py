import sys
class Cat:
    def __init__(self, files):
        self.files = files
        self.i = 0
        self.f = open(files[0], "rb")

    def read(self, n=-1):
        out = b""
        while n < 0 or len(out) < n:
            chunk = self.f.read(-1 if n < 0 else n - len(out))
            if chunk:
                out += chunk
                continue
            self.i += 1
            if self.i >= len(self.files):
                break
            self.f.close()
            self.f = open(self.files[self.i], "rb")
        return out


