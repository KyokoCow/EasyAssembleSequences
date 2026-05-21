import os

# =========================================================
# ログ出力関数のダミー。colabではこれ使う。
# =========================================================
def log_message(message):

    print(message)
# =========================================================
# Colab判定
# =========================================================

def is_colab():

    try:
        import google.colab
        return True

    except ImportError:
        return False


# =========================================================
# 各処理（後で中身を追加）
# =========================================================

def convert_seq_to_fasta(input_dir):

    converted_count = 0

    # .seq ファイル処理
    for filename in os.listdir(input_dir):

        if filename.endswith(".seq"):

            seq_path = os.path.join(
                input_dir,
                filename
            )

            # 拡張子除去
            sample_name = os.path.splitext(
                filename
            )[0]

            # 配列読み込み
            with open(
                seq_path,
                "r",
                encoding="utf-8"
            ) as f:

                sequence = f.read().strip()

            # FASTA形式
            fasta_text = (
                f">{sample_name}\n"
                f"{sequence}\n"
            )

            # 出力ファイル名
            fasta_filename = (
                sample_name + ".fasta"
            )

            # ★同じフォルダへ保存
            fasta_path = os.path.join(
                input_dir,
                fasta_filename
            )

            # 保存
            with open(
                fasta_path,
                "w",
                encoding="utf-8"
            ) as f:

                f.write(fasta_text)

            converted_count += 1

            # ログ
            log_message(
                f"変換完了: {fasta_filename}"
            )

    # 完了
    log_message(
        f"seq → fasta変換完了 "
        f"({converted_count} files)"
    )

def trim_fasta(input_dir):

    trimmed_count = 0

    # フォルダ内処理
    for filename in os.listdir(input_dir):

        if (
            filename.endswith(".seq")
            and (
                filename.startswith("F_")
                or filename.startswith("R_")
            )
        ):

            fasta_path = os.path.join(
                input_dir,
                filename
            )

            sample_name = os.path.splitext(
                filename
            )[0]

            # =====================================
            # FASTA読み込み
            # =====================================

            sequence_lines = []

            with open(
                fasta_path,
                "r",
                encoding="utf-8"
            ) as f:

                for line in f:

                    line = line.strip()

                    # ヘッダー除外
                    if line.startswith(">"):
                        continue

                    sequence_lines.append(line)

            sequence = "".join(sequence_lines)

            # =====================================
            # トリミング
            # =====================================

            trimmed_sequence = sequence[20:700]

            # =====================================
            # FASTA形式
            # =====================================

            fasta_text = (
                f">{sample_name}\n"
                f"{trimmed_sequence}\n"
            )

            # =====================================
            # 上書き保存
            # =====================================

            with open(
                fasta_path,
                "w",
                encoding="utf-8"
            ) as f:

                f.write(fasta_text)

            trimmed_count += 1

            # ログ
            log_message(
                f"トリミング完了: {filename}"
            )

    # =========================================
    # 完了
    # =========================================

    log_message(
        f"トリミング完了 "
        f"({trimmed_count} files)"
    )


def check_files(input_dir):

    # =====================================
    # ファイル一覧取得
    # =====================================

    files = sorted(
        os.listdir(input_dir)
    )

    # =====================================
    # ログ
    # =====================================

    log_message(
        f"ファイルチェック完了 "
        f"({len(files)} files)"
    )

    # =====================================
    # Colab表示
    # =====================================

    print("\n========== FILE LIST ==========")

    for file in files:

        print(file)

    # =====================================
    # GUIへ返す
    # =====================================

    return files
#######################################################################################
#ここからアセンブリ処理
#######################################################################################
print("ライブラリ読み込み中...")

import os
import re

# =========================
# 基本関数
# =========================

def read_sequence(filepath):

    seq = ""

    with open(filepath, "r", encoding="utf-8") as f:

        for line in f:

            line = line.strip().upper()

            if line.startswith(">") or line == "":
                continue

            # ATCGN以外除去
            line = re.sub(r"[^ATCGN]", "", line)

            seq += line

    return seq


def reverse_complement(seq):

    comp = str.maketrans(
        "ATCGN",
        "TAGCN"
    )

    return seq.translate(comp)[::-1]


def count_match(a, b):

    """
    -同士は一致にカウントしない
    """

    return sum(
        1
        for x, y in zip(a, b)
        if x == y and x != "-"
    )


def find_best_position_F(ref, query):

    """
    Forward用
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


def find_best_position_R(ref, query):

    """
    Reverse用
    """

    best_score = -1
    best_pos = 0

    q_len = len(query)
    r_len = len(ref)

    for i in range(0, r_len + q_len - 1):

        score = 0

        for j in range(q_len):

            ref_idx = i - (q_len - 1 - j)

            if 0 <= ref_idx < r_len:

                if ref[ref_idx] == query[j]:

                    score += 1

        if score > best_score:

            best_score = score
            best_pos = i

    return best_pos, best_score


# =========================
# アセンブリ本体
# =========================

def assembly_sequences(input_dir):

    data_dir = input_dir

    files = [

        os.path.join(data_dir, f)

        for f in os.listdir(data_dir)

        if f.endswith((".fasta", ".fa", ".seq"))
    ]

    ref_seq = None
    ref_name = "REF"

    forward_seq_data = []
    reverse_seq_data = []

    # =========================
    # ファイル分類
    # =========================

    for filepath in files:

        filename = os.path.basename(filepath)

        seq = read_sequence(filepath)

        # REF
        if filename.startswith("REF_"):

            ref_seq = seq

            # 拡張子除去
            ref_name = os.path.splitext(
                filename
            )[0]

        # Reverse
        elif filename.startswith("R_"):

            seq = reverse_complement(seq)

            reverse_seq_data.append(
                (filename, seq)
            )

        # Forward
        elif filename.startswith("F_"):

            forward_seq_data.append(
                (filename, seq)
            )

    # =========================
    # REF確認
    # =========================

    if ref_seq is None:

        log_message(
            "REF_ で始まる参照配列が見つかりません"
        )

        return

    log_message(
        f"参照配列長: {len(ref_seq)}"
    )

    # =========================
    # ソート
    # =========================

    forward_seq_data.sort(

        key=lambda x: int(
            re.search(
                r"F_(\d+)",
                x[0]
            ).group(1)
        )
    )

    reverse_seq_data.sort(

        key=lambda x: int(
            re.search(
                r"R_(\d+)",
                x[0]
            ).group(1)
        )
    )

    # =========================
    # Forward配置
    # =========================

    forward_placements = []

    for name, seq in forward_seq_data:

        pos, score = find_best_position_F(
            ref_seq,
            seq
        )

        log_message(
            f"[F] {name}: "
            f"position={pos}, "
            f"match={score}"
        )

        forward_placements.append(
            (name, seq, pos)
        )

    # =========================
    # Reverse配置
    # =========================

    reverse_placements = []

    for name, seq in reverse_seq_data:

        pos, score = find_best_position_R(
            ref_seq,
            seq
        )

        log_message(
            f"[R] {name}: "
            f"position={pos}, "
            f"match={score}"
        )

        reverse_placements.append(
            (name, seq, pos)
        )

    # =========================
    # Reverse補正
    # =========================

    reverse_placements_fixed = []

    for name, seq, pos in reverse_placements:

        start_pos = pos - (len(seq) - 1)

        reverse_placements_fixed.append(
            (name, seq, start_pos)
        )

    # =========================
    # 全配置統合
    # =========================

    all_placements = (

        forward_placements +
        reverse_placements_fixed
    )

    min_pos = min(

        [0] +

        [
            pos
            for _, _, pos
            in all_placements
        ]
    )

    max_end = max(

        [len(ref_seq)] +

        [
            pos + len(seq)
            for _, seq, pos
            in all_placements
        ]
    )

    total_len = max_end - min_pos

    # =========================
    # 出力
    # =========================

    output_file = os.path.join(
        data_dir,
        "aligned.fasta"
    )

    with open(
        output_file,
        "w",
        encoding="utf-8"
    ) as f:

        # =====================
        # REF
        # =====================

        aligned_ref = ["-"] * total_len

        for i, base in enumerate(ref_seq):

            idx = i - min_pos

            aligned_ref[idx] = base

        f.write(f">{ref_name}\n")

        for i in range(0, total_len, 60):

            f.write(
                "".join(
                    aligned_ref[i:i+60]
                ) + "\n"
            )

        # =====================
        # 全配列出力
        # =====================

        for name, seq, pos in all_placements:

            aligned_seq = ["-"] * total_len

            for i, base in enumerate(seq):

                idx = pos + i - min_pos

                if 0 <= idx < total_len:

                    aligned_seq[idx] = base

            header = os.path.splitext(name)[0]

            f.write(f">{header}\n")

            for i in range(
                0,
                total_len,
                60
            ):

                f.write(
                    "".join(
                        aligned_seq[i:i+60]
                    ) + "\n"
                )

    # =========================
    # 完了
    # =========================

    log_message(
        f"アセンブリ完了: {output_file}"
    )

#######################################################################################
#ここからコンセンサス処理
#######################################################################################
print("ライブラリ読み込み中...")

from collections import Counter
import os


# =========================
# FASTA読み込み
# =========================

def read_fasta(filepath):

    seqs = []

    name = None

    seq = []

    with open(
        filepath,
        "r",
        encoding="utf-8"
    ) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    seqs.append(
                        (
                            name,
                            "".join(seq)
                        )
                    )

                name = line[1:]

                seq = []

            else:

                seq.append(
                    line.upper()
                )

        if name is not None:

            seqs.append(
                (
                    name,
                    "".join(seq)
                )
            )

    return seqs


# =========================
# コンセンサス計算
# =========================

def build_consensus(aligned_seqs):

    ref_len = len(aligned_seqs[0][1])

    consensus = []

    confidence = []

    for i in range(ref_len):

        bases = []

        # "-" と "N" を除外
        for _, seq in aligned_seqs:

            b = seq[i]

            if b != "-" and b != "N":

                bases.append(b)

        # 全欠損
        if len(bases) == 0:

            consensus.append("N")

            confidence.append("T")

            continue

        counter = Counter(bases)

        base = counter.most_common(1)[0][0]

        consensus.append(base)

        # -------------------------
        # 一致率計算
        # -------------------------

        total = 0
        match = 0

        for _, seq in aligned_seqs:

            b = seq[i]

            if b == "-":
                continue

            total += 1

            if b == base:
                match += 1

        ratio = match / total if total > 0 else 0

        # -------------------------
        # confidence
        # -------------------------

        if total == 1:

            confidence.append("A")

        elif ratio == 1.0:

            confidence.append("C")

        elif ratio <= 0.5:

            confidence.append("T")

        else:

            confidence.append("G")

    return (
        "".join(consensus),
        "".join(confidence)
    )


# =========================
# FASTA保存
# =========================

def save_fasta(

    filepath,
    seqs,
    consensus,
    confidence
):

    with open(
        filepath,
        "w",
        encoding="utf-8"
    ) as f:

        # 元配列
        for name, seq in seqs:

            f.write(f">{name}\n")

            for i in range(
                0,
                len(seq),
                60
            ):

                f.write(
                    seq[i:i+60] + "\n"
                )

        # CONSENSUS
        f.write(">CONSENSUS\n")

        for i in range(
            0,
            len(consensus),
            60
        ):

            f.write(
                consensus[i:i+60] + "\n"
            )

        # CONFIDENCE
        f.write(">CONFIDENCE\n")

        for i in range(
            0,
            len(confidence),
            60
        ):

            f.write(
                confidence[i:i+60] + "\n"
            )


# =========================
# コンセンサス処理統合
# =========================

def consensus_pipeline(

    input_dir,
    mode="make"
):

    # =========================
    # モード設定
    # =========================

    if mode == "make":

        input_file = os.path.join(
            input_dir,
            "aligned.fasta"
        )

        output_file = os.path.join(
            input_dir,
            "aligned_with_consensus.fasta"
        )

    elif mode == "fix":

        input_file = os.path.join(
            input_dir,
            "aligned_with_consensus.fasta"
        )

        # 上書き保存
        output_file = input_file

    else:

        log_message(
            f"不明なmode: {mode}"
        )

        return

    # =========================
    # 入力ファイル確認
    # =========================

    if not os.path.exists(input_file):

        log_message(
            f"入力ファイルが見つかりません: "
            f"{input_file}"
        )

        return

    log_message(
        f"入力ファイル: {input_file}"
    )

    # =========================
    # FASTA読み込み
    # =========================

    aligned_seqs = read_fasta(
        input_file
    )

    # =========================
    # 除外対象
    # =========================

    aligned_seqs = [

        (name, seq)

        for name, seq in aligned_seqs

        if not (
        name.startswith("REF_")
        or name == "CONSENSUS"
        or name == "CONFIDENCE"
    )
    ]

    # =========================
    # コンセンサス計算
    # =========================

    consensus, confidence = build_consensus(
        aligned_seqs
    )

    # =========================
    # 保存
    # =========================

    save_fasta(

        output_file,
        aligned_seqs,
        consensus,
        confidence
    )

    # =========================
    # 完了
    # =========================

    log_message(
        f"コンセンサス処理完了: "
        f"{output_file}"
    )

# =========================================================
# Colab版UI
# =========================================================

def run_colab_ui():

    print("===================================")
    print(" Sequence Tool (Colab Mode)")
    print("===================================")

    input_dir = input(
        "\n対象フォルダを入力してください: "
    ).strip()

    if not input_dir:

        print("フォルダが入力されませんでした")

        return

    while True:

        print("\n-----------------------------------")
        print("実行する処理を選択してください")
        print("-----------------------------------")
    
        print("1 : フォルダ内ファイルチェック")
        print("2 : トリミング")
        print("3 : アセンブリ")
        print("4 : コンセンサス生成")
        print("5 : コンセンサス再生成")
        print("6 :（おまけ） seq → fasta変換")
        print("0 : 終了")

        choice = input(
            "\n番号入力: "
        ).strip()

        # =========================
        # seq → fasta
        # =========================

        if choice == "6":

            convert_seq_to_fasta(
                input_dir
            )

        # =========================
        # trimming
        # =========================

        elif choice == "2":

            trim_fasta(
                input_dir
            )

        # =========================
        # file check
        # =========================

        elif choice == "1":

            check_files(
                input_dir
            )

        # =========================
        # assembly
        # =========================

        elif choice == "3":

            assembly_sequences(
                input_dir
            )

        # =========================
        # consensus make
        # =========================

        elif choice == "4":

            consensus_pipeline(

                input_dir,

                mode="make"
            )

        # =========================
        # consensus fix
        # =========================

        elif choice == "5":

            consensus_pipeline(

                input_dir,

                mode="fix"
            )

        # =========================
        # end
        # =========================

        elif choice == "0":

            print("終了します")

            break

        # =========================
        # invalid
        # =========================

        else:

            print("無効な入力です")


# =========================================================
# ローカルGUI
# =========================================================

def run_local_gui():

    import tkinter as tk
    from tkinter import filedialog
    from tkinter import messagebox

    # --------------------------------
    # ログ出力関数を上書き
    # --------------------------------

    global log_message

    def log_message(message):

        log_text.insert(tk.END, message + "\n")

        log_text.see(tk.END)

    # --------------------------------
    # フォルダ参照
    # --------------------------------

    def browse_folder():

        folder = filedialog.askdirectory(
            title="データフォルダを選択"
        )

        if folder:
            folder_var.set(folder)

    # --------------------------------
    # 入力フォルダ取得
    # --------------------------------

    def get_input_dir():

        input_dir = folder_var.get().strip()

        if not input_dir:

            messagebox.showerror(
                "エラー",
                "データフォルダを選択してください"
            )

            return None

        return input_dir
    
    # =====================================================
    # ファイル一覧更新
    # =====================================================

    def update_file_list(files):

        # 一旦クリア
        file_listbox.delete(0, tk.END)

        # ファイル追加
        for file in files:

            file_listbox.insert(
                tk.END,
                file
            )

    # --------------------------------
    # 各ボタン処理
    # --------------------------------

    def run_convert():

        input_dir = get_input_dir()

        if not input_dir:
            return

        log_message("seq → fasta変換開始")

        convert_seq_to_fasta(input_dir)

        log_message("seq → fasta変換完了")


    def run_trim():

        input_dir = get_input_dir()

        if not input_dir:
            return

        log_message("トリミング開始")

        trim_fasta(input_dir)

        log_message("トリミング完了")


    def run_check():

        input_dir = get_input_dir()

        if not input_dir:
            return

        files = check_files(input_dir)

        update_file_list(files)


    def run_assembly():

        input_dir = get_input_dir()

        if not input_dir:
            return
        
        log_message("アセンブリ開始")

        assembly_sequences(input_dir)

        log_message("アセンブリ完了")


    def run_consensus():

        input_dir = get_input_dir()

        if not input_dir:
            return
        
        log_message("コンセンサス生成開始")

        consensus_pipeline(
            input_dir,
            mode="make"
        )

        log_message("コンセンサス生成完了")


    def run_fix():

        input_dir = get_input_dir()

        if not input_dir:
            return
        
        log_message("コンセンサス再生成開始")

        consensus_pipeline(
            input_dir,
            mode="fix"
        )

        log_message("コンセンサス再生成完了")

    
    # =====================================================
    # GUI作成
    # =====================================================

    root = tk.Tk()

    root.title("Sequence Tool")

    root.geometry("900x600")

    # フォルダ変数
    folder_var = tk.StringVar()

    # =====================================================
    # ログ出力関数
    # =====================================================

    def log_message(message):

        log_text.insert(tk.END, message + "\n")

        # 自動スクロール
        log_text.see(tk.END)

    # =====================================================
    # メインフレーム
    # =====================================================

    main_frame = tk.Frame(root)
    main_frame.pack(fill="both", expand=True)

    # -------------------------
    # 上部エリア
    # -------------------------

    top_frame = tk.Frame(main_frame)
    top_frame.pack(fill="both", expand=True)

    # -------------------------
    # 下部ログエリア
    # -------------------------

    log_frame = tk.Frame(main_frame)
    log_frame.pack(fill="x", padx=10, pady=10)

    # =====================================================
    # 左側：ファイル一覧
    # =====================================================

    left_frame = tk.Frame(top_frame)
    left_frame.pack(
        side="left",
        fill="both",
        expand=True,
        padx=10,
        pady=10
    )

    tk.Label(
        left_frame,
        text="ファイル一覧"
    ).pack(anchor="w")

    file_listbox = tk.Listbox(
        left_frame,
        width=50,
        height=25
    )

    file_listbox.pack(
        fill="both",
        expand=True
    )

    # スクロールバー
    scrollbar = tk.Scrollbar(file_listbox)

    scrollbar.pack(
        side="right",
        fill="y"
    )

    file_listbox.config(
        yscrollcommand=scrollbar.set
    )

    scrollbar.config(
        command=file_listbox.yview
    )

    # =====================================================
    # 右側：操作パネル
    # =====================================================

    right_frame = tk.Frame(top_frame)
    right_frame.pack(
        side="right",
        fill="y",
        padx=10,
        pady=10
    )

    # --------------------------------
    # フォルダ選択
    # --------------------------------

    tk.Label(
        right_frame,
        text="データフォルダ"
    ).pack(anchor="w")

    folder_entry = tk.Entry(
        right_frame,
        textvariable=folder_var,
        width=40
    )

    folder_entry.pack(pady=5)

    browse_button = tk.Button(
        right_frame,
        text="参照",
        width=20,
        command=browse_folder
    )

    browse_button.pack(pady=5)

    # --------------------------------
    # ボタン群
    # --------------------------------

    tk.Label(
        right_frame,
        text="\n実行する処理"
    ).pack()

    tk.Button(
        right_frame,
        text="フォルダ内ファイルチェック",
        width=30,
        command=run_check
    ).pack(pady=5)


    tk.Button(
        right_frame,
        text="トリミング",
        width=30,
        command=run_trim
    ).pack(pady=5)

    tk.Button(
        right_frame,
        text="アセンブリ",
        width=30,
        command=run_assembly
    ).pack(pady=5)

    tk.Button(
        right_frame,
        text="コンセンサス生成",
        width=30,
        command=run_consensus
    ).pack(pady=5)

    tk.Button(
        right_frame,
        text="コンセンサス再生成",
        width=30,
        command=run_fix
    ).pack(pady=5)

    tk.Button(
        right_frame,
        text="(おまけ) seq → fasta変換",
        width=30,
        command=run_convert
    ).pack(pady=5)

    # --------------------------------
    # 終了ボタン
    # --------------------------------

    tk.Button(
        right_frame,
        text="終了",
        width=20,
        command=root.destroy
    ).pack(pady=20)

    # =====================================================
    # ログ欄
    # =====================================================

    tk.Label(
        log_frame,
        text="ログ"
    ).pack(anchor="w")

    log_text = tk.Text(
        log_frame,
        height=8
    )

    log_text.pack(
        fill="x"
    )

    # 初期メッセージ
    log_message("Sequence Tool 起動")

    root.mainloop()


# =========================================================
# メイン
# =========================================================

if __name__ == "__main__":
        

        # Colab
        if is_colab():

            run_colab_ui()

        # ローカル
        else:

            run_local_gui()