# proyek1-daa-paleogenomika
# Proyek 1 Analisis Algoritma - Kasus A: Paleogenomika

Nama: Celsi Nabiva  |  NIM: 12550125685 |  Teknik Informatika, UIN Suska Riau

Program Python (OOP) yang memetakan fragmen DNA (16 basa) ke 100 sekuens spesies dinosaurus
(masing-masing 256 basa). Sebuah fragmen dianggap cocok jika identik 100%, atau beda maksimal
1 basa di ujung paling kanan (mutasi). Hasil diurutkan dari skor tertinggi.

## Algoritma
1. **KMP** mencari 15 karakter pertama pola secara eksak (tanpa `find()`, `index()`, `re`).
2. Karakter ke-16 dibandingkan: sama = *Exact Match* (16/16), beda = *End-Mutated* (15/16).
3. **Merge Sort** (rekursif dan iteratif) mengurutkan kandidat dan hasil akhir
   (skor menurun, indeks menaik). Tanpa `sorted()` / `.sort()`.
4. Setiap perbandingan karakter dihitung lewat kelas `Penghitung` sebagai metrik operasi.

## Struktur repositori
```
paleogenomika.py             kode utama (Spesies, PemuatDataset, Penghitung, KMP,
                             MergeSort, HasilPemetaan, PemetaDNA)
simulasi_paleogenomika.py    simulasi dan pengujian (terbaik, terburuk, skala input)
dataset/case1_dino_dna.csv   dataset Kasus A
README.md                    berkas ini
```

## Cara menjalankan
```
python paleogenomika.py dataset/case1_dino_dna.csv              # versi iteratif
python paleogenomika.py dataset/case1_dino_dna.csv --rekursif   # versi rekursif
python simulasi_paleogenomika.py dataset/case1_dino_dna.csv     # simulasi dan pengujian
```
Butuh Python 3.8 atau lebih baru; tidak memakai library eksternal.

## Contoh keluaran (3 pola contoh)
```
1. Pola: "CTTCGGGCAACTTAAC"  -> Plateosaurus trossingensis, nomor 45, indeks 32, 100.0% (Exact Match)
2. Pola: "CTGGCACAAACCGTTG"  -> Mapusaurus roseae, nomor 90, indeks 220, 93.8% (End-Mutated)
3. Pola: "TTTTTTTTTTTTTTTT"  -> tidak ditemukan (Not Found)
Total Operasi Dasar (Perbandingan Karakter): 95885 operasi
```
