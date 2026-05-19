# 参照配列ベース簡易アラインメントスクリプト

概要

本スクリプトは、複数のシーケンス断片を参照配列に基づいて**簡易的に位置合わせ（アラインメント）**するためのツールです。

各シーケンスは参照配列上でスライドさせ、一致塩基数が最大となる位置に配置されます。結果は、すべての配列を共通座標系上に並べたFASTA形式として出力されます（多重アラインメント風の表示）。

- フォワード／リバース配列に対応
- 参照配列の前後にはみ出す配列にも対応
- ギャップ（挿入・欠失）は考慮しない簡易アルゴリズム

---

# EasyAssembleSequences

サンガーシーケンスデータ用の簡易アセンブリ支援ツールです。

以下の環境に対応しています。

- ローカルGUI版（tkinter）
- Google Colab版

GitHub Repository:  
:contentReference[oaicite:0]{index=0}

---

# 主な機能

- seq → FASTA変換
- 配列トリミング
- フォルダ内ファイル確認
- 簡易アセンブリ
- コンセンサス配列生成
- コンセンサス再計算（修正用）

---

# ローカルでの使用方法

```bash
python easy_assemble.py
```

実行するとGUIが起動します。

---

# Google Colabでの使用方法

## リポジトリ取得

```python
!git clone https://github.com/KyokoCow/EasyAssembleSequences.git
```

## フォルダ移動

```python
%cd EasyAssembleSequences
```

## 実行

```python
!python easy_assemble.py
```

---

# Google Drive使用例

```python
from google.colab import drive
drive.mount('/content/drive')
```

例：

```text
/content/drive/MyDrive/data
```

を入力フォルダとして使用できます。

---

# 想定ファイル名

## 参照配列

```text
REF_*.fasta
```

## Forward配列

```text
F_*.seq
F_*.fasta
```

## Reverse配列

```text
R_*.seq
R_*.fasta
```

---

# 出力ファイル

| 処理 | 出力 |
|---|---|
| アセンブリ | aligned.fasta |
| コンセンサス生成 | aligned_with_consensus.fasta |
| 再計算 | RECHECK_aligned_with_consensus.fasta |

---

# 注意

簡易的な探索・確認用ツールとして作成しています。