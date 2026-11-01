"""
Simulasi dan pengujian Kasus A (Paleogenomika).

Menjalankan:
  1. Pengujian fungsional (pola contoh lembar) pada versi iteratif dan rekursif.
  2. Skenario terbaik dan terburuk (dataset sintetis 100 sekuens x 256 basa).
  3. Eksperimen skala input: operasi, waktu, dan memori puncak.

Cara menjalankan:
    python simulasi_paleogenomika.py case1_dino_dna.csv
"""

import random
import sys
import time
import tracemalloc

from paleogenomika import PemetaDNA, PemuatDataset, Spesies


def buat_dataset(jumlah, isi=None, panjang=256, seed=1):
    """Dataset sintetis. isi=None berarti acak A/C/G/T, selain itu 1 karakter berulang."""
    acak = random.Random(seed)
    daftar = []
    for i in range(1, jumlah + 1):
        if isi is None:
            seq = "".join(acak.choice("ACGT") for _ in range(panjang))
        else:
            seq = isi * panjang
        daftar.append(Spesies(i, f"Spesies {i}", seq))
    return daftar


def jalankan(dataset, patterns, rekursif):
    pemeta = PemetaDNA(dataset, rekursif=rekursif)
    awal = time.perf_counter()
    hasil = pemeta.petakan(patterns)
    detik = time.perf_counter() - awal
    return hasil, pemeta.penghitung.jumlah, detik


def memori_puncak(dataset, patterns, rekursif):
    tracemalloc.start()
    PemetaDNA(dataset, rekursif=rekursif).petakan(patterns)
    _, puncak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return puncak


def uji_fungsional(path):
    print("=== 1. PENGUJIAN FUNGSIONAL (pola contoh lembar) ===")
    dataset = PemuatDataset.muat(path)
    patterns = ["TTTTTTTTTTTTTTTT", "CTGGCACAAACCGTTG", "CTTCGGGCAACTTAAC"]
    harapan = [("Exact Match", 45, 32), ("End-Mutated", 90, 220), ("Not Found", "-", "-")]
    for rek in (False, True):
        hasil, ops, _ = jalankan(dataset, patterns, rek)
        sesuai = [(h.status, h.nomor, h.indeks) for h in hasil] == harapan
        print(f"{'Rekursif' if rek else 'Iteratif'}: sesuai harapan = {sesuai}, operasi = {ops}")


def uji_skenario():
    print("\n=== 2. SKENARIO TERBAIK DAN TERBURUK (100 sekuens x 256 basa) ===")
    skenario = [
        ("Terbaik: huruf pertama pola tidak pernah muncul",
         buat_dataset(100, isi="A"), ["C" * 16]),
        ("Terburuk A: banyak mundur (fallback) KMP",
         buat_dataset(100, isi="A"), ["A" * 14 + "C" + "G"]),
        ("Terburuk B: kandidat sangat banyak (semua posisi cocok)",
         buat_dataset(100, isi="A"), ["A" * 15 + "T"]),
    ]
    for nama, dataset, patterns in skenario:
        for rek in (False, True):
            hasil, ops, detik = jalankan(dataset, patterns, rek)
            mem = memori_puncak(dataset, patterns, rek)
            print(f"{nama}\n  {'Rekursif' if rek else 'Iteratif'}: status={hasil[0].status}, "
                  f"operasi={ops}, waktu={detik * 1000:.1f} ms, memori puncak={mem / 1024:.0f} KiB")


def uji_skala():
    print("\n=== 3. EKSPERIMEN SKALA INPUT (dataset acak, 3 pola) ===")
    patterns = ["TTTTTTTTTTTTTTTT", "CTGGCACAAACCGTTG", "CTTCGGGCAACTTAAC"]
    print("jumlah_sekuens | mode      | operasi | waktu (ms) | memori puncak (KiB)")
    for jumlah in (25, 50, 100, 200, 400):
        dataset = buat_dataset(jumlah)
        for rek in (False, True):
            _, ops, detik = jalankan(dataset, patterns, rek)
            mem = memori_puncak(dataset, patterns, rek)
            print(f"{jumlah:>14} | {'rekursif' if rek else 'iteratif':<9} | {ops:>7} | "
                  f"{detik * 1000:>10.1f} | {mem / 1024:>8.0f}")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "case1_dino_dna.csv"
    uji_fungsional(path)
    uji_skenario()
    uji_skala()
