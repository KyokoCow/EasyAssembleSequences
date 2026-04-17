import os

# =========================
# 基本関数
# =========================

def read_fasta(filepath):
    """FASTAを読み込んで配列（1本）を返す"""
    seq = ""
    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


def reverse_complement(seq):
    """逆相補鎖を返す"""
    comp = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(comp)[::-1]


def count_match(seq1, seq2):
    """一致塩基数をカウント（長さは同じ前提）"""
    return sum(1 for a, b in zip(seq1, seq2) if a == b)


def find_best_position(ref, query):
    """
    ref上でqueryをスライドさせて最も一致する位置を探す
    戻り値: (best_pos, best_score)
    """
    best_score = -1
    best_pos = 0

    q_len = len(query)
    r_len = len(ref)

    # スライド
    for i in range(r_len - q_len + 1):
        ref_sub = ref[i:i + q_len]
        score = count_match(ref_sub, query)

        if score > best_score:
            best_score = score
            best_pos = i

    return best_pos, best_score


# =========================
# メイン処理
# =========================

def main():
    files = [f for f in os.listdir() if f.endswith(".fasta") or f.endswith(".fa")]

    ref_seq = None
    seq_data = []

    # ファイル分類
    for f in files:
        if f.startswith("REF_"):
            ref_seq = read_fasta(f)
        else:
            seq = read_fasta(f)

            # R_なら逆相補
            if f.startswith("R_"):
                seq = reverse_complement(seq)

            seq_data.append((f, seq))

    if ref_seq is None:
        print("REF_ で始まる参照配列が見つかりません")
        return

    print("参照配列長:", len(ref_seq))

    # 各配列の最適位置を探索
    placements = []

    for name, seq in seq_data:
        pos, score = find_best_position(ref_seq, seq)
        print(f"{name}: position={pos}, match={score}")

        placements.append((name, seq, pos))

    # =========================
    # コンセンサス生成
    # =========================

    # 最大長を確保
    max_len = len(ref_seq)

    # 各位置の塩基を格納
    consensus = list(ref_seq)  # 初期は参照

    for name, seq, pos in placements:
        for i, base in enumerate(seq):
            ref_index = pos + i

            if ref_index >= max_len:
                continue

            # 参照がNなら上書き
            if consensus[ref_index] == "N":
                consensus[ref_index] = base

            # 一致ならそのまま
            elif consensus[ref_index] == base:
                continue

            # 不一致は単純にquery優先（簡易仕様）
            else:
                consensus[ref_index] = base

    consensus_seq = "".join(consensus)

    # =========================
    # 出力
    # =========================

    output_file = "assembled.fasta"

    with open(output_file, "w") as f:
        f.write(">assembled_sequence\n")

        # 60文字ごとに改行
        for i in range(0, len(consensus_seq), 60):
            f.write(consensus_seq[i:i+60] + "\n")

    print(f"\n出力完了: {output_file}")


# =========================
if __name__ == "__main__":
    main()