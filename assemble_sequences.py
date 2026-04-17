import os
import re

# =========================
# 基本関数
# =========================

def read_sequence(filepath):
    seq = ""
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip().upper()

            if line.startswith(">") or line == "":
                continue

            # ATCGN以外を除去
            line = re.sub(r"[^ATCGN]", "", line)

            seq += line

    return seq


def reverse_complement(seq):
    comp = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(comp)[::-1]


def count_match(a, b):
    return sum(1 for x, y in zip(a, b) if x == y)


def find_best_position(ref, query):
    """
    負の位置も含めて最適位置を探索
    """
    best_score = -1
    best_pos = 0

    q_len = len(query)
    r_len = len(ref)

    for i in range(-q_len + 1, r_len):
        score = 0

        for j in range(q_len):
            ref_idx = i + j

            if 0 <= ref_idx < r_len:
                if ref[ref_idx] == query[j]:
                    score += 1

        if score > best_score:
            best_score = score
            best_pos = i

    return best_pos, best_score


# =========================
# メイン処理
# =========================

def main():
    files = [f for f in os.listdir() if f.endswith((".fasta", ".fa", ".seq"))]

    ref_seq = None
    seq_data = []

    # ファイル分類
    for f in files:
        if f.startswith("REF_"):
            ref_seq = read_sequence(f)
        else:
            seq = read_sequence(f)

            if f.startswith("R_"):
                seq = reverse_complement(seq)

            seq_data.append((f, seq))

    if ref_seq is None:
        print("REF_ で始まる参照配列が見つかりません")
        return

    print("参照配列長:", len(ref_seq))

    # 配置計算
    placements = []

    for name, seq in seq_data:
        pos, score = find_best_position(ref_seq, seq)
        print(f"{name}: position={pos}, match={score}")
        placements.append((name, seq, pos))

    # =========================
    # 全体座標の決定（ここが重要）
    # =========================

    min_pos = min([0] + [pos for _, _, pos in placements])
    max_end = max([len(ref_seq)] + [pos + len(seq) for _, seq, pos in placements])

    total_len = max_end - min_pos

    # =========================
    # 出力
    # =========================

    output_file = "aligned.fasta"

    with open(output_file, "w") as f:

        # REF
        aligned_ref = ["-"] * total_len

        for i, base in enumerate(ref_seq):
            idx = i - min_pos
            aligned_ref[idx] = base

        f.write(">REF\n")
        for i in range(0, total_len, 60):
            f.write("".join(aligned_ref[i:i+60]) + "\n")

        # 各配列
        for name, seq, pos in placements:

            aligned_seq = ["-"] * total_len

            for i, base in enumerate(seq):
                idx = pos + i - min_pos
                if 0 <= idx < total_len:
                    aligned_seq[idx] = base

            f.write(f">{name}\n")
            for i in range(0, total_len, 60):
                f.write("".join(aligned_seq[i:i+60]) + "\n")

    print(f"\n出力完了: {output_file}")


# =========================
if __name__ == "__main__":
    main()