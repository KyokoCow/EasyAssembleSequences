参照配列ベース簡易アラインメントスクリプト

概要

本スクリプトは、複数のシーケンス断片を参照配列に基づいて**簡易的に位置合わせ（アラインメント）**するためのツールです。

各シーケンスは参照配列上でスライドさせ、一致塩基数が最大となる位置に配置されます。結果は、すべての配列を共通座標系上に並べたFASTA形式として出力されます（多重アラインメント風の表示）。

- フォワード／リバース配列に対応
- 参照配列の前後にはみ出す配列にも対応
- ギャップ（挿入・欠失）は考慮しない簡易アルゴリズム

---

使用方法

1. 以下のファイルを同一フォルダに配置します
   
   - 本スクリプト（"assemble_sequences.py"）
   - 参照配列（FASTA形式、ファイル名は "REF_" で始める）
   - シーケンスファイル（".fasta" / ".fa" / ".seq"）

2. コマンドプロンプトで実行：

python assemble_sequences.py

3. 出力ファイル：

aligned.fasta

---

ファイル命名ルール

- "REF_*.fasta"
  → 参照配列

- "R_*.fasta" または "R_*.seq"
  → 自動的に逆相補鎖に変換して処理

- その他
  → フォワード配列として処理

---

出力形式

出力ファイル（"aligned.fasta"）には：

- 参照配列（REF）
- 各シーケンス
- 共通座標系で配置された配列
- ギャップは "-" で表示

例：

>REF
ATGCTAGCTAGCTAG

>sample1
---CTAGCTAG-----

>sample2
CTAGCT----------

---

アルゴリズム

1. 各シーケンスについて：
   
   - 参照配列上でスライド
   - 負の位置（左側へのはみ出し）も含めて探索
   - 一致塩基数をカウント

2. 最も一致数の多い位置を採用

3. 全配列をカバーする共通座標系を構築：
   
   - 左右に必要に応じて拡張

4. 全配列をその座標系に配置して出力

---

特徴

- 参照配列の外側に伸びる配列にも対応
- "R_" プレフィックスで逆相補鎖を自動処理
- 外部ライブラリ不要
- 軽量・高速

---

制限事項

- ギャップ（挿入・欠失）を考慮しない
- 1ファイル1配列のみ対応
- スコアリングは単純な一致数のみ（ミスマッチペナルティなし）

---

動作環境

- Python 3.x
- 外部ライブラリ不要

---

想定用途

本スクリプトは、以下のような用途を想定しています：

- Sangerシーケンス断片の簡易比較
- 迅速な目視確認
- 軽量なスクリプトベースの解析

---

注意

本ツールは簡易的なアラインメント用です。
厳密な解析には、BLAST、MAFFT、MUSCLEなどの専用ツールの使用を推奨します。

Simple Sequence Alignment Script (Reference-based)

Overview

This script performs a simple, reference-based alignment of multiple sequence fragments.

Each sequence is aligned to a reference sequence by sliding it along the reference and selecting the position with the highest number of matching bases. The result is exported in a FASTA format where all sequences (including the reference) are displayed in a common coordinate system, similar to a multiple sequence alignment.

- Supports forward and reverse-complement sequences
- Handles sequences extending before or after the reference
- Does not consider gaps (insertions/deletions)

---

Usage

1. Place the following files in the same directory:
   
   - This script ("assemble_sequences.py")
   - Reference sequence (FASTA format, filename must start with "REF_")
   - Sequence files (".fasta", ".fa", or ".seq")

2. Run the script:

python assemble_sequences.py

3. Output file:

aligned.fasta

---

File Naming Rules

- "REF_*.fasta"
  → Reference sequence

- "R_*.fasta" or "R_*.seq"
  → Automatically converted to reverse-complement

- Others
  → Treated as forward sequences

---

Output Format

The output ("aligned.fasta") contains:

- The reference sequence
- All input sequences aligned to the reference
- A shared coordinate system
- Gaps represented as "-"

Example:

>REF
ATGCTAGCTAGCTAG

>sample1
---CTAGCTAG-----

>sample2
CTAGCT----------

---

Algorithm

1. For each sequence:
   
   - Slide across the reference sequence
   - Evaluate all positions including negative offsets (left extension)
   - Count matching bases

2. Select the position with the highest match score

3. Build a global coordinate system:
   
   - Extend left/right as needed to include all sequences

4. Output all sequences aligned to this coordinate

---

Features

- Supports sequences extending beyond the reference boundaries
- Automatic reverse-complement handling ("R_" prefix)
- No external dependencies
- Lightweight and fast

---

Limitations

- No gap-aware alignment (no insertion/deletion handling)
- Single sequence per file only
- Simple scoring (match count only, no mismatch penalty)

---

Requirements

- Python 3.x
- No external libraries required

---

Intended Use

This script is designed for quick and rough alignment of Sanger sequencing fragments, especially when:

- You need a fast visual comparison
- Full alignment tools are unnecessary
- You want a lightweight, script-based solution

---

Disclaimer

This tool is intended for approximate alignment only.
For rigorous sequence analysis, use dedicated alignment tools such as BLAST, MAFFT, or MUSCLE.