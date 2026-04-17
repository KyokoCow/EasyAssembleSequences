# シーケンス簡易アセンブリスクリプト

## 概要

複数のシーケンスデータを参照配列に基づいて位置合わせし、1本の配列として結合する簡易スクリプトです。
ギャップは考慮せず、一致塩基数が最大となる位置に各配列を配置します。

---

## 使用方法

1. 以下のファイルを同じフォルダに配置します

   * 本スクリプト（assemble_sequences.py）
   * 参照配列（FASTA形式、ファイル名は `REF_` で始める）
   * シーケンスデータ（FASTAまたは.seq形式）

2. コマンドプロンプトで以下を実行します

```
python assemble_sequences.py
```

3. 実行後、`assembled.fasta` が出力されます

---

## ファイル命名ルール

* `REF_***.fasta`
  → 参照配列

* `R_***.fasta` または `R_***.seq`
  → 逆相補鎖として処理される

* その他
  → 通常のシーケンスとして処理

---

## アルゴリズム

* 各シーケンスを参照配列上でスライドさせる
* 一致塩基数が最大となる位置を探索
* 最適位置に配置
* すべての配列を統合してコンセンサス配列を生成

---

## 注意事項

* ギャップ（挿入・欠失）は考慮していません
* 複数配列を含むFASTAには対応していません（1ファイル1配列）
* 不一致は単純に後から処理した配列を優先します

---

## 動作環境

* Python 3.x（外部ライブラリ不要）

---

## 補足

簡易的なアセンブリを目的としたスクリプトのため、
厳密な解析には専用ツールの使用を推奨します。


# Simple Sequence Assembly Script

## Overview

This script performs a simple assembly of multiple sequence fragments based on a reference sequence.
Each sequence is aligned by sliding along the reference and placed at the position with the highest number of matching bases.

This is a lightweight approach and does **not** consider gaps (insertions/deletions).

---

## Usage

1. Place the following files in the same directory:

   * This script (`assemble_sequences.py`)
   * Reference sequence (FASTA format, filename must start with `REF_`)
   * Sequence files (FASTA or `.seq` format)

2. Run the script from the command line:

```bash
python assemble_sequences.py
```

3. After execution, an output file named `assembled.fasta` will be generated.

---

## File Naming Rules

* `REF_*.fasta`
  → Reference sequence

* `R_*.fasta` or `R_*.seq`
  → Treated as reverse-complement sequences

* Others
  → Treated as forward sequences

---

## Algorithm

* Each sequence is slid along the reference sequence
* The position with the highest number of matching bases is selected
* Sequences are placed at their optimal positions
* A consensus sequence is generated from all placements

---

## Notes

* Insertions and deletions (gaps) are **not considered**
* Multi-sequence FASTA files are **not supported** (one sequence per file)
* In case of mismatches, later sequences overwrite earlier ones

---

## Requirements

* Python 3.x
* No external libraries required

---

## Disclaimer

This script is intended for **simple and quick assembly tasks only**.
For accurate sequence analysis, dedicated alignment tools are recommended.
