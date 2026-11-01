"""
Proyek 1 - Kasus A: Paleogenomika (Sekuens DNA Kuno)

Memetakan fragmen DNA (n = 16) ke 100 sekuens spesies (tiap 256 basa).
Pencocokan sah: identik 100% atau beda maksimal 1 karakter di ujung kanan.

Strategi:
  1. Cari 15 karakter pertama pola secara eksak dengan KMP.
  2. Bandingkan karakter ke-16: sama = Exact Match, beda = End-Mutated.
  3. Urutkan dengan Merge Sort (skor menurun, indeks menaik).

Tidak memakai str.find(), str.index(), re, sorted(), atau .sort().
Setiap perbandingan karakter dihitung lewat kelas Penghitung.

Cara menjalankan:
    python paleogenomika.py case1_dino_dna.csv
    python paleogenomika.py case1_dino_dna.csv --rekursif
"""

import csv
import sys


class Spesies:
    """Satu baris dataset: id, nama latin, dan sekuens DNA."""

    def __init__(self, id_spesies, nama_latin, sekuens):
        self.id = id_spesies
        self.nama_latin = nama_latin
        self.sekuens = sekuens


class PemuatDataset:
    """Membaca file CSV menjadi daftar objek Spesies."""

    @staticmethod
    def muat(path):
        daftar = []
        with open(path, newline="", encoding="utf-8") as f:
            for baris in csv.DictReader(f):
                daftar.append(
                    Spesies(int(baris["id"]), baris["nama_latin"], baris["sekuens_dna"])
                )
        return daftar


class Penghitung:
    """Membungkus perbandingan karakter dan menghitung total operasi dasar."""

    def __init__(self):
        self.jumlah = 0

    def bandingkan(self, a, b):
        self.jumlah += 1
        return a == b

    def reset(self):
        self.jumlah = 0


class KMP:
    """Algoritma Knuth-Morris-Pratt, versi iteratif dan rekursif."""

    def __init__(self, penghitung):
        self.penghitung = penghitung

    # ---------- versi non-rekursif ----------
    def bangun_failure(self, pola):
        m = len(pola)
        failure = [0] * m
        j = 0
        for i in range(1, m):
            while True:
                if self.penghitung.bandingkan(pola[i], pola[j]):
                    j += 1
                    break
                if j == 0:
                    break
                j = failure[j - 1]
            failure[i] = j
        return failure

    def cari(self, teks, pola, failure):
        m = len(pola)
        hasil = []
        j = 0
        for i in range(len(teks)):
            while True:
                if self.penghitung.bandingkan(teks[i], pola[j]):
                    j += 1
                    break
                if j == 0:
                    break
                j = failure[j - 1]
            if j == m:
                hasil.append(i - m + 1)
                j = failure[j - 1]
        return hasil

    # ---------- versi rekursif ----------
    def _lanjut(self, pola, failure, karakter, j):
        """Mencari nilai j baru setelah membaca satu karakter (rekursif)."""
        if self.penghitung.bandingkan(karakter, pola[j]):
            return j + 1
        if j == 0:
            return 0
        return self._lanjut(pola, failure, karakter, failure[j - 1])

    def bangun_failure_rekursif(self, pola):
        failure = [0] * len(pola)
        self._isi_failure(pola, failure, 1, 0)
        return failure

    def _isi_failure(self, pola, failure, i, j):
        if i == len(pola):
            return
        j = self._lanjut(pola, failure, pola[i], j)
        failure[i] = j
        self._isi_failure(pola, failure, i + 1, j)

    def cari_rekursif(self, teks, pola, failure):
        hasil = []
        self._pindai(teks, pola, failure, 0, 0, hasil)
        return hasil

    def _pindai(self, teks, pola, failure, i, j, hasil):
        if i == len(teks):
            return
        j = self._lanjut(pola, failure, teks[i], j)
        if j == len(pola):
            hasil.append(i - len(pola) + 1)
            j = failure[j - 1]
        self._pindai(teks, pola, failure, i + 1, j, hasil)


class MergeSort:
    """Merge Sort stabil, versi rekursif (top-down) dan iteratif (bottom-up)."""

    @staticmethod
    def _gabung(kiri, kanan, kunci):
        hasil = []
        i = j = 0
        while i < len(kiri) and j < len(kanan):
            if kunci(kiri[i]) <= kunci(kanan[j]):
                hasil.append(kiri[i])
                i += 1
            else:
                hasil.append(kanan[j])
                j += 1
        while i < len(kiri):
            hasil.append(kiri[i])
            i += 1
        while j < len(kanan):
            hasil.append(kanan[j])
            j += 1
        return hasil

    @staticmethod
    def urutkan_rekursif(daftar, kunci):
        if len(daftar) <= 1:
            return list(daftar)
        tengah = len(daftar) // 2
        kiri = MergeSort.urutkan_rekursif(daftar[:tengah], kunci)
        kanan = MergeSort.urutkan_rekursif(daftar[tengah:], kunci)
        return MergeSort._gabung(kiri, kanan, kunci)

    @staticmethod
    def urutkan_iteratif(daftar, kunci):
        blok = [[x] for x in daftar]
        while len(blok) > 1:
            baru = []
            for i in range(0, len(blok), 2):
                if i + 1 < len(blok):
                    baru.append(MergeSort._gabung(blok[i], blok[i + 1], kunci))
                else:
                    baru.append(blok[i])
            blok = baru
        return blok[0] if blok else []


class HasilPemetaan:
    """Satu baris laporan untuk satu pola."""

    def __init__(self, pola, nama_latin="-", nomor="-", indeks="-",
                 cocok=0, status="Not Found"):
        self.pola = pola
        self.nama_latin = nama_latin
        self.nomor = nomor
        self.indeks = indeks
        self.cocok = cocok
        self.status = status

    @property
    def skor(self):
        return 100.0 * self.cocok / len(self.pola)


class PemetaDNA:
    """Mengatur alur: cari awalan, cek karakter ujung, pilih terbaik, urutkan."""

    def __init__(self, dataset, rekursif=False):
        self.dataset = dataset
        self.rekursif = rekursif
        self.penghitung = Penghitung()
        self.kmp = KMP(self.penghitung)

    def _petakan_satu(self, pola):
        n = len(pola)
        awalan = pola[:-1]
        if self.rekursif:
            failure = self.kmp.bangun_failure_rekursif(awalan)
        else:
            failure = self.kmp.bangun_failure(awalan)

        kandidat = []
        for spesies in self.dataset:
            if self.rekursif:
                posisi = self.kmp.cari_rekursif(spesies.sekuens, awalan, failure)
            else:
                posisi = self.kmp.cari(spesies.sekuens, awalan, failure)
            for s in posisi:
                if s + n > len(spesies.sekuens):
                    continue
                sama = self.penghitung.bandingkan(spesies.sekuens[s + n - 1], pola[-1])
                kandidat.append(HasilPemetaan(
                    pola, spesies.nama_latin, spesies.id, s,
                    n if sama else n - 1,
                    "Exact Match" if sama else "End-Mutated",
                ))

        if not kandidat:
            return HasilPemetaan(pola)
        urut = self._urutkan(kandidat)
        return urut[0]

    @staticmethod
    def _kunci(h):
        indeks = h.indeks if h.indeks != "-" else 10 ** 9
        nomor = h.nomor if h.nomor != "-" else 10 ** 9
        return (-h.cocok, indeks, nomor)

    def _urutkan(self, daftar):
        if self.rekursif:
            return MergeSort.urutkan_rekursif(daftar, self._kunci)
        return MergeSort.urutkan_iteratif(daftar, self._kunci)

    def petakan(self, patterns):
        self.penghitung.reset()
        hasil = [self._petakan_satu(p) for p in patterns]
        return self._urutkan(hasil)

    def cetak_laporan(self, hasil):
        print("=== HASIL PEMETAAN FRAGMEN DNA (TERURUT SKOR TERTINGGI) ===")
        for no, h in enumerate(hasil, start=1):
            n = len(h.pola)
            skor = "0%" if h.cocok == 0 else f"{h.skor:.1f}%"
            print(f'\n{no}. Pola: "{h.pola}"')
            print(f"   Nama Latin: {h.nama_latin}")
            print(f"   Nomor: {h.nomor}")
            print(f"   Indeks: {h.indeks}")
            print(f"   Skor Kemiripan: {skor} ({h.cocok}/{n})")
            print(f"   Status: {h.status}")
        print("\nMetrik Eksekusi:")
        print(f"- Total Operasi Dasar (Perbandingan Karakter): {self.penghitung.jumlah} operasi")


def main():
    path = "case1_dino_dna.csv"
    rekursif = False
    for arg in sys.argv[1:]:
        if arg == "--rekursif":
            rekursif = True
        else:
            path = arg

    patterns = [
        "TTTTTTTTTTTTTTTT",
        "CTGGCACAAACCGTTG",
        "CTTCGGGCAACTTAAC",
    ]
    dataset = PemuatDataset.muat(path)
    pemeta = PemetaDNA(dataset, rekursif=rekursif)
    hasil = pemeta.petakan(patterns)
    pemeta.cetak_laporan(hasil)


if __name__ == "__main__":
    main()
